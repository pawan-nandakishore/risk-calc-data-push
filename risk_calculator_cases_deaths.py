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
print(file_name)

# Initalize the s3 client to push data to s3 
s3_client = boto3.client('s3', aws_access_key_id=access_token, aws_secret_access_key=secret_key)


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

def smooth_data(group):
        filled_cases = group['ConfirmedCases'].fillna(0)
        smoothed_cases = gaussian_filter1d(filled_cases.values,7)
        group['Smooth7ConfirmedCases'] = smoothed_cases
        
        filled_deaths = group['ConfirmedDeaths'].fillna(0)
        smoothed_deaths = gaussian_filter1d(filled_deaths.values,7)
        group['Smooth7ConfirmedDeaths'] = smoothed_deaths
        return group
    
def get_cases_deaths(country_subdivisons): 
    data_loc = "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
    oxford_df = pd.read_csv(data_loc)
    us_states = country_subdivisons['US']
    us_states = [x.replace("-", "_" ) for x in us_states]
    us_states_only = oxford_df[oxford_df['RegionCode'].isin(us_states)]
    states_subset = us_states_only[['Date','RegionCode', 'RegionName',  'ConfirmedCases', 'ConfirmedDeaths']]
    smoothed_states = states_subset.groupby("RegionCode").apply(lambda x : smooth_data(x))
    return smoothed_states
    




if __name__ == "__main__": 
    
    start_time = time.time()
    # Get population data
    population_df = pd.read_csv(file_name)

    # Get dictionary of states for a given country
    country_subdivisons = get_country_states()

    # Get today's date
    today = datetime.today().strftime('%Y-%m-%d')

    # Get smoothed confirmed cases and deaths data for all US states
    all_data = get_cases_deaths(country_subdivisons)

    # Save path for the S3 file 
    save_full_path = 'processed/risk-calculator-data/USStates_confirmed_cases_deaths.csv'

    # This function pushes the data to S3. Refer to the docstring for more details
    df_to_s3(all_data, s3_client, bucket_name, save_full_path)

    end_time = time.time() - start_time 

    print("Total number of seconds {}".format(end_time))
