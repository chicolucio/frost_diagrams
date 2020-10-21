from frost import potentials
import pandas as pd
from pandas.testing import assert_frame_equal
from ast import literal_eval


DATA_FILE = '../data/frost.csv'


def database(path):
    try:
        df = pd.read_csv(path)
        df = df.set_index('charge')
        return df
    except IOError:
        print('Cannot open data file')


DF = database(DATA_FILE)


def test_Cr_standard_acid():
    element = 'Cr'
    temperature = 298.15
    conc_ion = 1
    pH = 0
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_acid_bal'
    df_golden = database('test_data_Cr_01.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)


def test_Cr_acid_change_concentration():
    element = 'Cr'
    temperature = 298.15
    conc_ion = 0.0005
    pH = 0
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_acid_bal'
    df_golden = database('test_data_Cr_02.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)


def test_Cr_acid_change_pH():
    element = 'Cr'
    temperature = 298.15
    conc_ion = 1
    pH = 3
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_acid_bal'
    df_golden = database('test_data_Cr_03.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)


def test_Cr_acid_change_temperature():
    element = 'Cr'
    temperature = 310
    conc_ion = 1
    pH = 0
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_acid_bal'
    df_golden = database('test_data_Cr_04.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)


def test_Cr_standard_basic():
    element = 'Cr'
    temperature = 298.15
    conc_ion = 1
    pH = 14
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_basic_bal'
    df_golden = database('test_data_Cr_05.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)


def test_Cr_basic_change_concentration():
    element = 'Cr'
    temperature = 298.15
    conc_ion = 0.0005
    pH = 14
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_basic_bal'
    df_golden = database('test_data_Cr_06.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)


def test_Cr_basic_change_pH_temperature():
    element = 'Cr'
    temperature = 310
    conc_ion = 1
    pH = 10
    df_potentials = potentials(element, DF, temperature=temperature, conc_ion=conc_ion, pH=pH)
    column_coefficients = element + '_basic_bal'
    df_golden = database('test_data_Cr_07.csv')
    df_golden[column_coefficients] = df_golden[column_coefficients].apply(literal_eval)
    return assert_frame_equal(df_golden, df_potentials)
