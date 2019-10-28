#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:32:24 2019

@author: nicholas
"""

import numpy as np


# An implementation of rao's spacing test but instead of using degrees as the
# unit I use the parameter N to set where wraparound happens. See
# http://jammalam.faculty.pstat.ucsb.edu/html/favorite/test.htm for a
# description of the method
def rao_spacing(data, N):
    m = data.size
    sorted_data = np.sort(data)
    diff = np.zeros_like(sorted_data)
    diff[:-1] = np.diff(sorted_data)
    diff[-1] = N + sorted_data[0] - sorted_data[-1]
    deviation_from_avg_spacing = np.abs(diff - N/m)
    return (180/N)*np.sum(deviation_from_avg_spacing)


def censored_sample(data, intervals):
    '''Return an array of the elements of data that do not fall within
    the intervals'''
    cens_sample = data
    for interval in intervals:
        to_remove = np.where(np.logical_and(data >= interval[0],
                                            data <= interval[1]))
        cens_sample = np.delete(cens_sample, to_remove)
    return cens_sample


# taken from stack_overflow question
# https://stackoverflow.com/questions/15273693/python-union-of-multiple-ranges
def union_intervals(intervals):
    b = []
    for begin, end in sorted(intervals):
        if b and (b[-1][1] >= begin - 1):
            b[-1] = (b[-1][0], end)
        else:
            b.append((begin, end))
    return b
