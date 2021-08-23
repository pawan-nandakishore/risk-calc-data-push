import pandas as pd
from datetime import datetime
import requests
import json 
import pycountry as pc
import time 
import os 
from scipy.ndimage import gaussian_filter1d
from dotenv import load_dotenv
import boto3 
import io
load_dotenv()
access_token = os.getenv("covid_s3_access_key")
secret_key = os.getenv("covid_s3_token")
bucket_name = os.getenv("bucket_name")
file_name = os.getenv("data_file")


# Main source for the training data
DATA_URL = 'https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv'
# Latest Oxford data
DATA_FILE = 'processed/risk-calculator-data/OxCGRT_latest.csv'


# Initalize the s3 client to push data to s3 
s3_client = boto3.client('s3', aws_access_key_id=access_token, aws_secret_access_key=secret_key)

# Get the oxford file
obj = s3_client.get_object(Bucket=bucket_name, Key=DATA_FILE)    
oxford_data = pd.read_csv(io.BytesIO(obj['Body'].read()), 
                    parse_dates=['Date'],
                    encoding="ISO-8859-1",
                    dtype={"RegionName": str,
                            "RegionCode": str},
                    error_bad_lines=False)

#### FUNCTIONS  ##########
def df_to_s3(df, s3_client, bucket_name, save_full_path):
    """ Covert a dataframe to csv and send to s3 

    Parameters
    ----------
    df : Pandas dataframe
        The input pandas dataframe to be sent.
    s3_client : boto3_client
        AWS boto3 client with s3 as the service
        where credentials have already been added.
    bucket_name : str
        S3 bucket name
    save_full_path : str
        The path in S3 for saving the file. 
    """

    # Only run this if the dataframe is not empty
    # else print error message and return nothing
    if not df.empty:                                     

        with io.StringIO() as csv_buffer:                          
            df.to_csv(csv_buffer, index=False)          # Convert dataframe to csv using csv buffer 
            
            # put object in S3, key is the name of file in s3
            # body is the csv buffer 
            response = s3_client.put_object(
                                            Bucket=bucket_name,
                                            Key=save_full_path,
                                            Body=csv_buffer.getvalue()
                                        )

            # Get status of the put operation. If status is 200 means the 
            # operation was successful else you print its unsuccessful with 
            # status code                            
            status = response.get("ResponseMetadata", {}).get("HTTPStatusCode")
            if status == 200:
                print("Successful S3 put_object response. Status - {}".format(status))
            else:
                print("Unsuccessful S3 put_object response. Status - {}".format(status))
    else: 
        print("Data frame does not exist hence no push was made")
    return 



def get_country_states():
    """ Generate a dict that countries the subdivisions 
        of a country and store them in a dictionary.
    Returns
    -------
    dict 
        Keys are countries and values are lists.
        Each list contains the subdivision code for the country.
    """
    all_subdivisions = list(pc.subdivisions)                    # Convert all subdivisions to list
    all_countries = [x.alpha_2 for x in list(pc.countries)]     # Store the alpha_2 code for each country 

    # Initialize a dictionary with empty list which we will populate with states
    countries_subdivisons = {key: [] for key in all_countries}

    # Add subdivision code to countries_subdivision based on the country code 
    for subdivision in all_subdivisions: 
            countries_subdivisons[subdivision.country_code].append(subdivision.code)
    
    return countries_subdivisons


def get_strains_world(start_date = '2020-01-01', end_day = None):
    
    new_strains_df = pd.DataFrame()

    dates = pd.date_range(start=start_date,end=end_day)
    dates = pd.DataFrame({'Date':dates})
    
    new_strains_df = pd.DataFrame()
    n_days = len(dates)

    #Get lineages WORLD
    for location in oxford_data.CountryCode.unique():
        country_name = oxford_data.loc[oxford_data.CountryCode == location, 'CountryName']
        response = requests.get(f'https://api.outbreak.info/genomics/prevalence-by-location-all-lineages?location_id={location}&ndays={n_days}').text
        if(json.loads(response)["success"]):
            results =  json.loads(response)['results']
            results_df = pd.DataFrame.from_dict(results)
            location_strains_df = dates.copy()
            population = population_df.loc[population_df.CountryCode == location,'Population'].max()
            location_strains_df['CountryCode'] = location
            location_strains_df['CountryName'] = country_name.iloc[0]
            location_strains_df['Population'] = population
            for lineage in results_df.lineage.unique():
                df_lineage = results_df.loc[results_df.lineage == lineage]
                new_strain = pd.DataFrame()
                new_strain['Date'] = pd.to_datetime(df_lineage['date'])
                prevalence = df_lineage['prevalence_rolling'].fillna(0)
                location_strains_df  = location_strains_df.merge(new_strain, how='left')
                prevalence14 = prevalence.fillna(0).rolling(14).mean()
                location_strains_df[f'prevalence_gaussian5_{lineage}'] = pd.Series(gaussian_filter1d(prevalence14.fillna(0),5))
            new_strains_df = new_strains_df.append(location_strains_df, ignore_index=True)
            print("Reading data for {}".format(location))
        else:
            print(f'¸No data available for {location}')
            
    return new_strains_df



def get_states_strains(country, country_subdivisons, start_date = '2020-01-01', end_day = None):
    

    alpha2_code = pc.countries.lookup(country).alpha_2
    state_codes = country_subdivisons[alpha2_code]
    alpha3_code = pc.countries.lookup(country).alpha_3

    new_strains_df = pd.DataFrame()

    dates = pd.date_range(start=start_date,end=end_day)
    dates = pd.DataFrame({'Date':dates})
    
    n_days = len(dates)

    #Get lineages WORLD
    for location in list(state_codes):
        state_name = pc.subdivisions.lookup(location).name 
        response = requests.get(f'https://api.outbreak.info/genomics/prevalence-by-location-all-lineages?location_id={alpha3_code}_{location}&ndays={n_days}').text
        if(json.loads(response)["success"]):
            print(f'Reading data for {location}')
            results =  json.loads(response)['results']
            results_df = pd.DataFrame.from_dict(results)
            location_strains_df = dates.copy()
            location_strains_df['CountryCode'] = location
            location_strains_df['CountryName'] = state_name
            for lineage in results_df.lineage.unique():
                df_lineage = results_df.loc[results_df.lineage == lineage]
                new_strain = pd.DataFrame()
                new_strain['Date'] = pd.to_datetime(df_lineage['date'])
                prevalence14 = df_lineage['prevalence_rolling'].fillna(0).rolling(14).mean()
                df_lineage[f'prevalence_gaussian5_{lineage}'] = gaussian_filter1d(prevalence14.fillna(0),5)
                df_lineage_subset = df_lineage[['date', f'prevalence_gaussian5_{lineage}']]
                df_lineage_subset['Date'] = df_lineage_subset['date'].apply(pd.to_datetime)
                df_lineage_subset = df_lineage_subset.drop(columns=["date"])
                location_strains_df  = pd.merge(location_strains_df, new_strain, left_on="Date", right_on="Date", how='left')
                location_strains_df = pd.merge(location_strains_df, df_lineage_subset, left_on="Date", right_on="Date", how='left')
            new_strains_df = new_strains_df.append(location_strains_df, ignore_index=True)
        else:
            print(f'!!! No data available for {location} !!!')
            
    return new_strains_df

if __name__ == "__main__": 
    
    start_time = time.time()
    # Get population data
    population_df = pd.read_csv(file_name)

    # Get dictionary of states for a given country
    country_subdivisons = get_country_states()


    # Get today's date
    today = datetime.today().strftime('%Y-%m-%d')

    # Pull data for all countries in the world
    world_data = get_strains_world(start_date = '2020-01-01', end_day=today)

    # Pull data for states of the countries in the countries list
    countries_list = ["United States of America", "India", "Brazil", "Canada"]
    data_dfs = []
    for country in countries_list:
        df = get_states_strains(country, country_subdivisons, start_date='2020-01-01', end_day=today)
        data_dfs.append(df)

    # Concat the states data
    all_data = pd.concat(data_dfs)

    # Concat the states data with the world data 
    all_data = pd.concat([world_data, all_data])

    # Save path for the S3 file 
    save_full_path = 'processed/risk-calculator-data/OxCGRT_latest.csv'

    # This function pushes the data to S3. Refer to the docstring for more details
    df_to_s3(all_data, s3_client, bucket_name, save_full_path)


    end_time = time.time() - start_time 

    print("Total number of seconds {}".format(end_time))