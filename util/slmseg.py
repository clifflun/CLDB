###############################################################################
#   SLMSuite - Python wrapper
#   Copyright (C) 2016  Valerio Orlandini, Alberto Magi
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


# Requires rpy2 package and Python 3
import rpy2.robjects as robjects
import sys
import csv

# Simplification of vector and matrix handling
from rpy2.robjects import FloatVector as Vector


def Matrix(values, rows):
    return robjects.r['matrix'](values, nrow=rows)

# SLMSeg R functions import and Python wrapping
from rpy2.robjects.packages import importr
SLMSeg = importr('SLMSeg')
Rbase = importr('base')

signal_ = list()
pos_ = list()
fw_ = int()
mw_ = 1
eta_ = float()
stepeta_ = int()
omega_ = float()
mi_ = float()
smu_ = float()
sepsilon_ = float()
muk_ = Vector([0.0])
total_pred_break_ = Vector([0.0])
total_pred_break_filtered_ = Vector([0.0])
data_seg_ = list()
        
def reset_seg_param():
    global signal_ 
    global pos_ 
    global fw_ 
    global mw_ 
    global eta_ 
    global stepeta_ 
    global omega_ 
    global mi_ 
    global smu_ 
    global sepsilon_ 
    global muk_ 
    global total_pred_break_ 
    global total_pred_break_filtered_ 
    global data_seg_ 

    signal_ = list()
    pos_ = list()
    fw_ = int()
    mw_ = 1
    eta_ = float()
    stepeta_ = int()
    omega_ = float()
    mi_ = float()
    smu_ = float()
    sepsilon_ = float()
    muk_ = Vector([0.0])
    total_pred_break_ = Vector([0.0])
    total_pred_break_filtered_ = Vector([0.0])
    data_seg_ = list()
    
def load_signal_file(filename, chromosome):
    infile = open(filename, 'r', newline='')
    intable = csv.reader(infile, delimiter='\t')
    global signal_
    global pos_
    for row in intable:
        chr = row[0]
        position = row[1]
        log2r = row[2]
        if chr.isdigit():
            if int(chr) == chromosome:
                signal_.append(float(log2r))
                pos_.append(int(position))
    signal_ = Matrix(Vector((signal_)), 1)
    pos_ = Vector((pos_))


def load_data(signal, pos):
    global signal_
    global pos_
    for i in signal:
        signal_.append(i)
    signal_ = Matrix(Vector((signal_)), 1)
    if len(signal) == len(pos):
        for i in pos:
            pos_.append(i)
        pos_ = Vector((pos_))
    else:
        print("Positions vector must have the same size of the Signals one")

def set_variables(omega, eta, stepeta, fw):
    global omega_
    global eta_
    global stepeta_
    global fw_
    omega_ = omega
    eta_ = eta
    stepeta_ = stepeta
    fw_ = fw


def param_est_seq():
    result = SLMSeg.ParamEstSeq(signal_, omega_)
    global mi_
    global smu_
    global sepsilon_
    mi_ = result[0]
    smu_ = result[1]
    sepsilon_ = result[2]


def muk_est():
    global muk_
    muk_ = SLMSeg.MukEst(signal_, mw_)


def joint_seg():
    global total_pred_break_
    total_pred_break_ = SLMSeg.JointSeg(signal_, eta_, omega_, muk_, mi_, smu_,
                                        sepsilon_)


def joint_seg_in():
    global total_pred_break_
    total_pred_break_ = SLMSeg.JointSegIn(signal_, eta_, omega_, muk_, mi_,
                                          smu_, sepsilon_, pos_, stepeta_)


def filter_seg():
    global total_pred_break_filtered_
    total_pred_break_filtered_ = SLMSeg.FilterSeg(total_pred_break_, fw_)


def seg_results():
    global data_seg_
    result = SLMSeg.SegResults(signal_, total_pred_break_filtered_)
    reslen = Rbase.length(result)[0]
    for i in range(0, reslen - 1):
        data_seg_.append(result[i])


def SLM():
    param_est_seq()
    muk_est()
    joint_seg()
    filter_seg()
    seg_results()


def HSLM():
    param_est_seq()
    muk_est()
    joint_seg_in()
    filter_seg()
    seg_results()


def signal():
    return signal_


def pos():
    return pos_


def fw():
    return fw_


def mw():
    return mw_


def eta():
    return eta_


def stepeta():
    return stepeta_


def omega():
    return omega_


def mi():
    return mi_


def smu():
    return smu_


def sepsilon():
    return sepsilon_


def muk():
    return muk_


def total_pred_break():
    return total_pred_break_


def total_pred_break_filtered():
    return total_pred_break_filtered_


def data_seg():
    return data_seg_
