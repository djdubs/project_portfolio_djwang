## Running R code

R code can be run to replicate the results in the [projects/results]{.underline} directory. Libraries may need to be installed prior to running the script. The required list of libraries is displayed at the beginning of each script.\
If a script uses external data, you will need to create a sub-folder [data]{.underline} under the [projects]{.underline} directory. The directions to obtain the data for each project are shown in the section below. For privacy reasons, the data for some projects may be excluded.

## Datasets

**Spatial Analysis**\
Link: https://data.wake.gov/search\
Note: To access datasets, search for

- Raleigh Police Incidents (NIBRS)
  - It will help to filter out all rows where City of Incident and District are blank or "NaN".
  - Blanks, missing values, or zeroes for latitude and longitude should also be filtered out.
- Townships
  - When downloading, choose the Shapefile option
