
from meteostat import Point, Daily
from datetime import datetime

# Define a list of locations with their corresponding dates (latitude, longitude, date)
locations_with_dates = [
(30.436734,114.193377,'2014-06-10','2024-04-05'),
]

for lat, lon, date_str1, date_str2 in locations_with_dates:
    # Parse the date strings into datetime objects
    date1 = datetime.strptime(date_str1, '%Y-%m-%d')
    date2 = datetime.strptime(date_str2, '%Y-%m-%d')
    
    # Define the location with altitude set to 0 (or any default value)
    location = Point(lat, lon, 0)
    
    # Loop through the two dates
    for date in [date1, date2]:
        # Define the start and end date as the same date to get data for a single day
        start = date
        end = date

        # Get daily data for the defined period
        data = Daily(location, start, end)
        data = data.fetch()

        # Print the data
        print(f"Location: ({lat}, {lon}), Date: {date.date()}")
        print(data)
