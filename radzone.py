#!/usr/bin/env python

"""
radzone.py: Plot radzone intervals over the HRC AntiCo Shield rate across
            the Chandra mission lifetime.
"""

__author__ = "Dr. Grant R. Tremblay"
__license__ = "MIT"

import os
import sys

import time
import datetime as dt

from astropy.io import ascii
from astropy.table import Table
from astropy.table import vstack

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num

import numpy as np
from scipy import stats


def parse_shieldrate(msid_directory, shieldrate_msid):

    # Make sure the .csv file exists before trying this:
    if os.path.isfile(msid_directory + shieldrate_msid + ".csv"):
        msid = ascii.read(msid_directory +
                          shieldrate_msid + ".csv",
                          format="fast_csv")

        print("Shieldrate MSID Parsed")
    else:
        print("MSID CSV file not present")
        sys.exit(1)

    rate = msid['midvals']
    time = convert_chandra_time(msid['times'])

    shieldrate = {"Time": time,
                  "Rate": rate}

    return shieldrate


def convert_chandra_time(rawtimes):
    """
    Convert input CXC time (sec) to the time base required for the matplotlib
    plot_date function (days since start of the Year 1 A.D - yes, really).
    :param times: iterable list of times, in units of CXCsec (sec since 1998.0)
    :rtype: plot_date times (days since Year 1 A.D.)
    """

    # rawtimes is in units of CXC seconds, or seconds since 1998.0
    # Compute the Delta T between 1998.0 (CXC's Epoch) and 1970.0 (Unix Epoch)

    seconds_since_1998_0 = rawtimes[0]

    cxctime = dt.datetime(1998, 1, 1, 0, 0, 0)
    unixtime = dt.datetime(1970, 1, 1, 0, 0, 0)

    # Calculate the first offset from 1970.0, needed by matplotlib's plotdate
    # The below is equivalent (within a few tens of seconds) to the command
    # t0 = Chandra.Time.DateTime(times[0]).unix
    delta_time = (cxctime - unixtime).total_seconds() + seconds_since_1998_0

    plotdate_start = epoch2num(delta_time)

    # Now we use a relative offset from plotdate_start
    # the number 86,400 below is the number of seconds in a UTC day

    chandratime = (np.asarray(rawtimes) -
                   rawtimes[0]) / 86400. + plotdate_start

    return chandratime


def styleplots():
    plt.style.use('ggplot')

    labelsizes = 13

    plt.rcParams['font.size'] = labelsizes
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = labelsizes
    plt.rcParams['xtick.labelsize'] = labelsizes
    plt.rcParams['ytick.labelsize'] = labelsizes


def makeplot(shieldrate):

    styleplots()

    fig, ax = plt.subplots(figsize=(12, 8))

    plt.axhline(y=65000, color='gray', alpha=0.8,
                label="SCS 107 Threshold (65,000 cps)")

    ax.plot_date(shieldrate["Time"], shieldrate["Rate"], markersize=0.8,
                 label='HRC Shield Rate (2SHEV1RT)')

    ax.set_yscale('log')
    ax.set_ylim(1000, 200000)
    ax.set_ylabel(r'Counts s$^{-1}$')
    ax.set_xlabel('Date')

    ax.legend()

    print("Plot constructed, showing now")

    plt.show()


def main():

    msid_directory = "msids/"
    shieldrate_msid = "2SHEV1RT"

    spacecraft_event_directory = "spacecraft_events/"
    spacecraft_event = "radzone_intervals"

    # This is a dictionary:
    shieldrate = parse_shieldrate(msid_directory, shieldrate_msid)

    makeplot(shieldrate)

if __name__ == '__main__':
    start_time = time.time()

    main()

    runtime = round((time.time() - start_time), 3)
    print("Finished in {} seconds".format(runtime))
