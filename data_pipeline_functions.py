import pandas as pd
import json 
from urllib.request import urlopen, HTTPError
import pycountry as ctr
import datetime
import os 
import boto3
import io 
from dotenv import load_dotenv
from scipy.ndimage import gaussian_filter1d
import time 
from tqdm import tqdm

load_dotenv()

access_token = os.getenv("covid_s3_access_key")
secret_key = os.getenv("covid_s3_token")
bucket_name = os.getenv("bucket_name")

s3_client = boto3.client(
    "s3",
    aws_access_key_id=access_token,
    aws_secret_access_key=secret_key,
   )




def get_json(source):
    """ Generate a json file a url.

    Parameters
    ----------
    source : str
        URL of the the file to be read.
    Returns
    -------
        json
        A json file generated using the json model.
    """
    try:
        handler = urlopen(source)
        status  = 0 
    except HTTPError as e:  
        handler = e.read()
        status = 1
        print("Failed to retrive data from URL")

    if status == 0: 
        with urlopen(source) as url:
            data = json.loads(url.read().decode()) 
        return data
    else: 
        return None

def get_country_states():
    """ Generate a dict that countries the subdivisions 
        of a country and store them in a dictionary.
    Returns
    -------
    dict 
        Keys are countries and values are lists.
        Each list contains the subdivision code for the country.
    """
    all_subdivisions = list(ctr.subdivisions)                    # Convert all subdivisions to list
    all_countries = [x.alpha_2 for x in list(ctr.countries)]     # Store the alpha_2 code for each country 

    # Initialize a dictionary with empty list which we will populate with states
    countries_subdivisons = {key: [] for key in all_countries}

    # Add subdivision code to countries_subdivision based on the country code 
    for subdivision in all_subdivisions: 
            countries_subdivisons[subdivision.country_code].append(subdivision.code)
    
    return countries_subdivisons


    
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
    
    

def generate_save_name(lineage, alpha3_code, folder_path="processed/variants", country=False):
    """ Generate a string which contains the save paths for saving S3 files

    Parameters
    ----------
    lineage : str
         Covid lineage for which folder is created
    alpha3_code : str
        Contains the alpha3_code for the country
    folder_path : str, optional
        [description], by default "processed/variants"
        folder structure for S3 
    Returns
    -------
        str, str
        Provide the full path and the save path. Full path is to the file 
        to be saved and save path is for the folder where file is saved
    """

    date_today = str(datetime.date.today())
    save_path = "{}/{}/{}/{}".format(folder_path, date_today, lineage, alpha3_code)
    if country == False:
        save_full_path = "{}/{}_lineage_data.csv".format(save_path, lineage)
    else:
        save_full_path = "{}/{}_lineage_data_country.csv".format(save_path, lineage)

    return save_full_path, save_path


def save_local(df, save_path, save_full_path):
    """ Save a dataframe locally

    Parameters
    ----------
    df : dataframe
        dataframe file to be saved
    save_path : str
        Folder path to be saved 
    save_full_path : str
        File path to be saved 
    """
    # Make directory is diretory does not exist
    if not os.path.exists(save_path):
            os.makedirs(save_path)
    df.to_csv(save_full_path, index=False)         # Save dataframe to csv file, do not save index
    return 

def states_smooth(df, column, sigma):
    """ Apply gaussian smoothing to a column in a
        dataframe. This function creates a column 
        after smooth. It does not modify the existing
        column 
    Parameters
    ----------
    df : dataframe
        Dataframe which contains column to be smoothed
    column: str
        Name of the column which 
    sigma : int
        Smoothing factor 

    Returns
    -------
    dataframe
        return dataframe after smooth operation has been performed
    """
    column_cleaned = "{}_cleaned".format(column)
    df[column_cleaned] = df[column].fillna(0)                                   # If there ar NaN values replace them with zeros
    df = df.sort_values(by="date")                                              # Sort values by date
    prevalence_col = "Smooth{}".format(sigma)                                   # Column name for smoothed data column 
    df[prevalence_col] = gaussian_filter1d(df[column_cleaned].values, sigma)    # Smooth the data using scipy function 
    return df

def pull_states_data(states_abbr, alpha3_code, lineage, states_dfs, no_data_states, sigma=7 ):
    """ Pull covid variant data for states of a country from outbreak.info

    Parameters
    ----------
    states_abbr : list
        List of states in a country
    alpha3_code : str
        Alpha3_code for a country
    lineage : str
        Covid variant 
    states_dfs : dict
        Dictonary to store dataframes for each state 
    no_data_states : list
        List to store any states for which there is 
        no variant data 
    sigma: int 

    Returns
    -------
    dict, list
        Output a dictionary of dataframes where each 
        dataframe represents a state. The second output
        is a list of states for which no data was available
    """

    # Loop through each state, pull data from outbreak.info, 
    # if the data exists then smooth the data
    for abbr in states_abbr: 
        pango_lineage =  "https://api.outbreak.info/genomics/prevalence-by-location?location_id={}_{}&pangolin_lineage={}".format(alpha3_code ,abbr, lineage)
        pango_data = get_json(pango_lineage)
        pango_states_df = pd.DataFrame(data=pango_data)['results'].apply(pd.Series)
        
        # If the dataframe exist then create a new column called abbr 
        # which contains the abbreviation for the state, the rename the 
        # proportion column as prevalence, then smooth the data and store
        # it in the dict states_dfs
        if pango_states_df.shape[0] != 0:
            pango_states_df['abbr'] = abbr
            pango_states_df = pango_states_df.rename(columns={"proportion":"prevalence"})
            pango_states_df = states_smooth(pango_states_df,column="prevalence", sigma=7)
            states_dfs[abbr]= pango_states_df
        else: 
            no_data_states.append(abbr)
    return states_dfs, no_data_states
    

def states_by_lineage(lineage, country, countries_states, s3_client, bucket_name, folder_path="processed/variants", export=False, local=False): 
    """ Pull COVID variant data for a given country and push it to S3

    Parameters
    ----------
    lineage : str
        Covid lineage
    country : str
        Normal name of the country, non-iso or alpha codes,
        check official spelling 
    countries_states : dict
        Dictionary containing the states for a all countries
    s3_client : aws_client
        AWS boto3 client with s3 as the service
        where credentials have already been added   
    bucket_name : str
        Name of S3 bucket where data is to be stored
    export : bool, optional
        Option to export the data, by default False
    local : bool, optional
        Option to export data, if False then data is sent
        to S3, if True then a local file is created
        ,by default False

    Returns
    -------
    dataframe, list
        If code is successful then returns a dataframe contains
        the covid data for all states of a country
    """

    # Initalizating variables 
    alpha3_code = ctr.countries.lookup(country).alpha_3
    alpha2_code = ctr.countries.lookup(country).alpha_2
    states_abbr = countries_states[alpha2_code] 
    states_dfs = {}
    no_data_states = []
    no_data_lineage_states = {}

    # Given the states of a country and the covid lineage pull the time series of covid 
    # cases for a given lineage
    states_dfs, no_data_states = pull_states_data(states_abbr, alpha3_code, lineage, states_dfs, no_data_states )
    
    # If the data exists then merged all the dataframes for all states into a single
    # dataframe and save it
    if list(states_dfs.values()):   
        bfts_states = pd.concat(list(states_dfs.values()))
        
        if export == True: 
            save_full_path, save_path = generate_save_name(lineage, alpha3_code, folder_path)               # Generate a save name for dataframe 
            
            if local==True:
                save_local(bfts_states, save_path, save_full_path)                             # Save dataframe locally
            else:   
                df_to_s3(bfts_states, s3_client, bucket_name, save_full_path)                  # Send data to S3
                no_data_lineage_states[lineage] = no_data_states
        return bfts_states, no_data_lineage_states
    
    else: 
        print("No data to save")
        return None, None 



def vaccination_daily_push(s3_client, bucket_name, folder_path="raw/vaccinations/daily"): 
    """ Get vaccination data from OWID github repository and push tables to S3.
        This data is to be pushed daily
    Parameters
    ----------
    s3_client : AWS_client
        AWS boto3 client with s3 as the service
        where credentials have already been added
    bucket_name : str
        Name of S3 bucket where data is stored
    folder_path: str 
        Folder where each file is to be stored
        by default folder_path= "raw/vaccinations/daily"
    """

    # URLs of tables which need to be pushed to S3
    locations = {"countries": "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv",
                 "age_group": "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations-by-age-group.csv",
                "manufacturers": "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations-by-manufacturer.csv",
                "us_states": "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv"}
    
    date_today = str(datetime.date.today())                                                 # Today's data for creating date based folder
    for name in locations.keys():   
        df = pd.read_csv(locations[name])                                                   # Create dataframe from US
        save_full_path = "{}/{}/{}".format(folder_path, date_today, name)                   # Specify file save name
        df_to_s3(df, s3_client, bucket_name, save_full_path )                               # Push file to S3
         
    return 


def vaccination_weekly_push(s3_client, bucket_name, folder_path="raw/vaccinations/weekly"):
    """ Get vaccination data from OWID. These tables are to be pushed weekly 
        
    Parameters
    ----------
    s3_client : AWS_client
        AWS boto3 client with s3 as the service
        where credentials have already been added
    bucket_name : str
        Name of S3 bucket where data is stored
    folder_path: str 
        Folder where each file is to be stored
        by default folder_path= "raw/vaccinations/weekly"
    """

    # operations as 
    locations = {"locations": "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/locations.csv"}
    date_today = str(datetime.date.today())
    for name in locations.keys(): 
        df = pd.read_csv(locations[name])
        save_full_path = "{}/{}/{}".format(folder_path, date_today, name)
        df_to_s3(df, s3_client, bucket_name, save_full_path )
    return 


def daily_data_cases_deaths(group): 
    """ Generate a new column that contain daily covid cases 
        and replace NaN values with zeros

    Parameters
    ----------
    group : dataframe
        A group of data from a pandas dataframe. 
        group is generated by the groupby operation
    Returns
    -------
    dataframe
        Modifying and returning the same group
    """

    # Take a diff of the confirmed cases and generate 
    # daily cases. Create a new column for daily cases
    # and store it. Fill missing values with 0 
    group['DailyCases'] = group['ConfirmedCases'].diff()
    group['DailyCases'] = group['DailyCases'].fillna(0)
    
    # Do the same as above but for deaths 
    group['DailyDeaths'] = group['ConfirmedDeaths'].diff()
    group['DailyDeaths'] = group['DailyDeaths'].fillna(0)
    return group

def smooth_cases_deaths(group, sigma):
    """ Apply smoothing function to a group of data from a dataframe

    Parameters
    ----------
    group : dataframe
         A group of data from a pandas dataframe. 
        group is generated by the groupby operation
    sigma : int
        Standard deviation for gaussian smoothing

    Returns
    -------
    dataframe
        Modifying and returning the same group

    """
    group = daily_data_cases_deaths(group)
    cases_col = "SmoothDailyCases{}".format(sigma)
    deaths_col = "SmoothDailyDeaths{}".format(sigma)
    
    group[cases_col] = gaussian_filter1d(group['DailyCases'].values, sigma)
    group[deaths_col] = gaussian_filter1d(group['DailyDeaths'].values, sigma)    
    return group



def oxford_smooth_and_push(s3_client, bucket_name, folder_path="processed/oxford_all"):
    """ Get the oxford dataset, smooth the daily cases and deaths 
        and push the data to S3 

    Parameters
    ----------
    s3_client : AWS_client
        AWS boto3 client with s3 as the service
        where credentials have already been added
    bucket_name : str
        Name of S3 bucket where data is stored
    folder_path : str, optional
        Folder in S3 where the files will be stored.
        By default "processed/oxford_all"
    """
    sigma = 7 
    date_today = str(datetime.date.today())
   
    # Get oxford data 
    all_data = pd.read_csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv")

    # Break the dataset into two parts National and state totals
    national = all_data[all_data['Jurisdiction']== "NAT_TOTAL"]
    states = all_data[all_data['Jurisdiction']== "STATE_TOTAL"]
    
    # Apply smoothing function to each group
    grouped_states = states.groupby(by=['CountryCode','RegionCode' ]).apply(lambda x: smooth_cases_deaths(x, sigma))
    grouped_national = national.groupby(by=['CountryCode']).apply(lambda x: smooth_cases_deaths(x, sigma))
    
    # Send files to S3
    files = {'states':grouped_states, 'national':grouped_national}
    for f in files.keys(): 
        save_full_path = "{}/{}/{}".format(folder_path, date_today,f)
        df_to_s3(files[f], s3_client, bucket_name, save_full_path)

    return

def country_by_lineage(lineage, country, s3_client, bucket_name, folder_path= "interim/variants", export=False, local=False): 
    """ Get data for a given country based on lineage, here we are not getting the data 
        by subdivision/state of the country     

    Parameters
    ----------
    lineage : str
        specify covid lineage in captial form  
    country : iso2 code of the country
        specify the iso2 code of the country
    s3_client : AWS_client
        S3 aws client
    bucket_name : str
        Name of the S3 bucket
    export : bool, optional
        Specify if you want to export the file or on, by default False
    local : bool, optional
        Specify where to export file, locally or onto S3. If set to False file will be exported to S3, by default False

    Returns
    -------
        dataframe 
        Returns the dataframe containing the country data
    """
    

    alpha3_code = ctr.countries.lookup(country).alpha_3
        
    pango_lineage =  "https://api.outbreak.info/genomics/prevalence-by-location?location_id={}&pangolin_lineage={}".format(alpha3_code, lineage)
    pango_data = get_json(pango_lineage)
    country_df = pd.DataFrame(data=pango_data)['results'].apply(pd.Series)
    country_df = country_df.rename(columns={"proportion":"prevalence"})
    country_df = states_smooth(country_df,column="prevalence", sigma=7)

    if list(country_df):   
        
        if export == True: 
            save_full_path, save_path = generate_save_name(lineage, alpha3_code, folder_path, country=True) 
                
            if local==True:
                save_local(country_df, save_path, save_full_path)
            else:
                df_to_s3(country_df,s3_client, bucket_name, save_full_path)
        return country_df
    else: 
        print("No data to save")
        return None






def s3_pull_latest_date(path, bucket): 
    """ Pull a file specified by path from a S3 bucket based 
    on the latest date. To do this we start with the latest 
    date try to get the data if it fails then go to the previous
    day and pull data if that fails then go to 2 days before
    and pull data, so on. 

    Parameters
    ----------
    path : str
        location of the file
    bucket : str
        Name of the S3 bucket

    Returns
    -------
    dict, int
        Outputs a dictionary that countains the location of individual files
    """
    
    num_days = 0  # This starts with Today's date
    keycount = 0  # Keycount is 0 if there is no data found
    
    # Try to pull data for the latest day, if it fails then try for yesterday
    # if that fails try for the day before yesterday, so on. It tries this for 
    # 10 days. 
    while keycount==0 or num_days >= 10:
        latest_date = str(datetime.date.today()-datetime.timedelta(days=num_days))
        response = s3_client.list_objects_v2(
                Bucket=bucket,
                Prefix ='{}/{}'.format(path,latest_date),
                MaxKeys=1000)
        keycount = response['KeyCount']
        num_days += 1

    # If the data is not available for the last 10 days 
    # then set response to none and keycount to 0 
    if num_days >= 10: 
        response = None 
        keycount = 0
        print("No data found for last 100 days. Please check code")
     
    return response, keycount

def  download_variants_data(bucket):
    """ Download variants data from S3 bucket

    Parameters
    ----------
    bucket : str
        S3 Bucket name
    """
    
    # Given the path of the variants folder pull data from 
    # S3 bucket
    path = "interim/variants"
    response, _ = s3_pull_latest_date(path, bucket)

    # If local folder with the name of path does not 
    # exist then create it    
    if not os.path.exists(path):
            os.makedirs(path)
    
    # Collect the paths of all files in the variants folder
    all_files = [key["Key"] for key in response['Contents']]
   
   # Pull file from S3 and save file into local folder
    for file in all_files:
        save_name = "_".join(file.split("/")[-2:])
        s3_client.download_file(bucket_name, file, "{}/{}".format(path,save_name))
    
    return 



def download_vaccine_data(bucket):
    """ Download vaccine data

    Parameters
    ----------
    bucket : str
        S3 Bucket name 
    """
    # Paths for daily and weekly data stores for the 
    # Oxford dataset
    save_paths = ["daily", "weekly"]
    for s in save_paths:     
        save_path ="vaccinations/{}".format(s)
        
        # If local folder with the name of path does not 
        # exist then create it   
        if not os.path.exists(save_path):
            os.makedirs(save_path) 

        # Pull paths to all files in the relevant folders in S3
        response, _ = s3_pull_latest_date("raw/{}".format(save_path), bucket)

        # Collect the paths of all files in the vaccinations folder
        all_files = [key["Key"] for key in response['Contents']]

        # Pull file from S3 and save file into local folder
        for file in all_files:
                save_name = "_".join(file.split("/")[-2:])
                s3_client.download_file(bucket_name, file, "{}/{}".format(save_path, save_name))
    
    return 

        

def download_oxford_data(bucket): 
    """ Download oxford dataset

    Parameters
    ----------
    bucket : str
        S3 bucket name
    """
    # Paths for daily and weekly data stores for the 
    # Oxford dataset
    path = "interim/oxford_all"

    # Pull paths to all files in the relevant folders in S3
    response, _ = s3_pull_latest_date(path, bucket)

    # If local folder with the name of path does not 
    # exist then create it   
    if not os.path.exists(path):
            os.makedirs(path)

    # Collect the paths of all files in the oxford_all folder
    all_files = [key["Key"] for key in response['Contents']]


    # Pull file from S3 and save file into local folder
    for file in all_files:
        save_name = "_".join(file.split("/")[-2:])
        s3_client.download_file(bucket_name, file, "{}/{}".format(path,save_name))
    return


def get_lineages(country, ndays):
    """ Get all lineages present in the country, based on the number of days

    Parameters
    ----------
    country : str
        Iso2 code of the country 
    ndays: int
        Number of days from current day to get the data variants 
    Returns
    -------
    list
        list of variants for the given country
    """
    alpha3_code = ctr.countries.lookup(country).alpha_3
    pango_lineage = 'https://api.outbreak.info/genomics/prevalence-by-location-all-lineages?location_id={}&ndays={}'.format(alpha3_code, ndays)
    pango_data = get_json(pango_lineage)
    if pango_data: 
        pango_states_df = pd.DataFrame(data=pango_data)['results'].apply(pd.Series)
        all_lineages = all_lineages = pango_states_df['lineage'].value_counts().index.tolist()
        if 'other' in all_lineages: 
            all_lineages.remove('other')

        all_lineages = [x.upper() for x in all_lineages]
    else: 
        all_lineages = []
    return all_lineages




if __name__ == "__main__": 

    start_time = time.time()

    countries_states = get_country_states()
    input_countries = list(countries_states.keys())
    ndays = 365
    for country in tqdm(input_countries): 
            
            alpha3_code = ctr.countries.lookup(country).alpha_3
            all_lineages = get_lineages(alpha3_code, ndays)
            
            for lineage in all_lineages:
                print("Country, lineage: {}, {}".format(alpha3_code, lineage))
                if  len(countries_states[country]) > 0 : 
                    _ , _ = states_by_lineage(lineage, 
                                                alpha3_code,
                                                countries_states,
                                                s3_client,
                                                bucket_name,
                                                folder_path="interim/variants",
                                                export=True,
                                                local=False)
        
            
                _= country_by_lineage(lineage, alpha3_code, s3_client, bucket_name, folder_path="interim/variants", export=True, local=False) 
        
                
    time_delta = time.time() - start_time
    print("Total time: {}".format(time_delta))


