from datetime import datetime as dt
import time

def toDtObject(date_string):
    # anticipates 'yyyy/mm/dd'
    return dt(int(date_string.split('/')[0]),
              int(date_string.split('/')[1]),
              int(date_string.split('/')[2]))

def toYearFraction_help(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch
    
    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)
    
    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration
    
    return date.year + fraction

def toYearFraction(date):
    return toYearFraction_help(toDtObject(date))
