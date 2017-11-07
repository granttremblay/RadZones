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


def parse_shieldrate(msid_directory, shieldrate_msid_filenames):

    bothfiles_found = os.path.isfile(msid_directory + shieldrate_msid_filenames["5min"]) and os.path.isfile(
        msid_directory + shieldrate_msid_filenames["Full"])

    # Make sure the .csv file exists before trying this:
    if bothfiles_found:
        fivemin = ascii.read(msid_directory +
                             shieldrate_msid_filenames["5min"],
                             format="fast_csv")

        full = ascii.read(msid_directory +
                          shieldrate_msid_filenames["Full"],
                          format="fast_csv")

        print("Shieldrate MSID parsed")
    else:
        print("I can't find both of the MSIDs needed")
        sys.exit(1)

    rate = fivemin['midvals']
    fullrate = full['vals']
    time = convert_chandra_time(fivemin['times'])
    fulltimes = convert_chandra_time(full['times'])

    shieldrate = {"Time": time,
                  "Rate": rate,
                  "Full Time": fulltimes,
                  "Full Rate": fullrate}

    return shieldrate


def parse_orbits(spacecraft_event_directory, spacecraft_event_filename):

    # Make sure the .csv file exists before trying this:
    if os.path.isfile(spacecraft_event_directory + spacecraft_event_filename):
        msid = ascii.read(spacecraft_event_directory +
                          spacecraft_event_filename)

        print("Spacecraft orbits parsed")
    else:
        print("MSID CSV file not present")
        sys.exit(1)

    # Available fields in Orbit table:
    # start,stop,tstart,tstop,dur,orbit_num,perigee,apogee,t_perigee,
    # start_radzone,stop_radzone,dt_start_radzone,dt_stop_radzone

    # Times are given like: 2000:003:15:27:47.271, so you need to convert
    # them into an mpl date.

    radzone_entry = convert_orbit_time(msid['start_radzone'])
    radzone_exit = convert_orbit_time(msid['stop_radzone'])

    orbit = {"Radzone Entry": radzone_entry,
             "Radzone Exit": radzone_exit}

    return orbit


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


def convert_orbit_time(rawtimes):
    """
    The orbit table gives times in the format: 2000:003:15:27:47.271, i.e.
    YYYY:DOY:HH:MM:SS.sss, so you need to convert these into a matplotlib date.
    """

    # Using %S.%f at the end converts to microseconds. I tested this
    # and it's correct.

    orbittime = []

    for i in range(len(rawtimes)):
        orbittime.append(dt.datetime.strptime(
            rawtimes[i], "%Y:%j:%H:%M:%S.%f"))

    return orbittime


def styleplots():
    plt.style.use('ggplot')

    labelsizes = 13

    plt.rcParams['font.size'] = labelsizes
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = labelsizes
    plt.rcParams['xtick.labelsize'] = labelsizes
    plt.rcParams['ytick.labelsize'] = labelsizes


def makeplot(shieldrate, orbit):

    styleplots()

    fig, ax = plt.subplots(figsize=(12, 8))

    # Fill in a gray bar for every radzone passage
    # and mark the entries and exits with a separate color bar

    # Yes, this is the dumb way you have to select ggplot colors:
    entry_color = list(plt.rcParams['axes.prop_cycle'])[2]['color']
    exit_color = list(plt.rcParams['axes.prop_cycle'])[4]['color']

    for i, (entry, exit) in enumerate(zip(orbit["Radzone Entry"], orbit["Radzone Exit"])):
        plt.axvline(entry, label="Radzone Entry" if i ==
                    0 else "", color=entry_color)
        plt.axvline(exit, label="Radzone Exit" if i ==
                    0 else "", color=exit_color)
        plt.axvspan(entry, exit, alpha=0.5, color='gray',
                    label="Radzone Passage" if i == 0 else "")
    # the label= business is to ensure I don't write thousands of entries in the legend

    # Also mark the entries and exits, to make it clear

    plt.axhline(y=65000, color='gray', alpha=0.8,
                label="SCS 107 Threshold (65,000 cps)")

    ax.plot_date(shieldrate["Time"], shieldrate["Rate"], markersize=0.8,
                 label='HRC Shield Rate (2SHEV1RT), 5 minute')

    ax.plot_date(shieldrate["Full Time"], shieldrate["Full Rate"], markersize=0.8,
                 label='HRC Shield Rate (2SHEV1RT), Full resolution')

    ax.set_yscale('log')
    ax.set_ylim(1000, 200000)
    ax.set_ylabel(r'Counts s$^{-1}$')
    ax.set_xlabel('Date')

    ax.legend()

    print("Plot constructed, showing now")

    plt.show()


def main():

    msid_directory = "msids/"
    #shieldrate_msid_filename = "2SHEV1RT_full.csv"
    shieldrate_msid_filenames = {"5min": "2SHEV1RT_5min.csv",
                                 "Full": "2SHEV1RT_full.csv"}

    spacecraft_event_directory = "spacecraft_events/"
    spacecraft_event_filename = "orbit_table.csv"

    # This is a dictionary:
    shieldrate = parse_shieldrate(msid_directory, shieldrate_msid_filenames)
    orbit = parse_orbits(spacecraft_event_directory, spacecraft_event_filename)

    print("There have been {} radzone passages since the start of the Kadi database".format(
        len(orbit["Radzone Entry"])))
    makeplot(shieldrate, orbit)


if __name__ == '__main__':
    start_time = time.time()

    main()

    runtime = round((time.time() - start_time), 3)
    print("Finished in {} seconds".format(runtime))
