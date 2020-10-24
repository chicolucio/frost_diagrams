def columns_df(elem_symbol):
    column_names = {'acid': ''.join((elem_symbol, '_acid')),
                    'basic': ''.join((elem_symbol, '_basic')),
                    'acid_state': ''.join((elem_symbol, '_acid', '_state')),
                    'basic_state': ''.join((elem_symbol, '_basic', '_state')),
                    'acid_label': ''.join((elem_symbol, '_acid_label')),
                    'basic_label': ''.join((elem_symbol, '_basic_label')),
                    'acid_bal': ''.join((elem_symbol, '_acid_bal')),
                    'basic_bal': ''.join((elem_symbol, '_basic_bal')),
                    'acid_E0': ''.join((elem_symbol, '_acid', '_E0')),
                    'basic_E0': ''.join((elem_symbol, '_basic', '_E0')),
                    'acid_E': ''.join((elem_symbol, '_acid', '_E')),
                    'basic_E': ''.join((elem_symbol, '_basic', '_E')),
                    'acid_new': ''.join((elem_symbol, '_acid', '_new')),
                    'basic_new': ''.join((elem_symbol, '_basic', '_new')),
                    }
    return column_names


def coef_ac(iterable):
    a, b, c, d, m = iterable
    x = 1
    z = c
    w = b / d
    y = 2 * z - a
    n = m + y
    return (x, y, n, w, z)


def coef_bas(iterable):
    a, b, c, d, m = iterable
    x = 1
    w = b / d
    z = c
    r = 2 * z - a
    s = r
    n = m + s
    return (x, r, n, w, z, s)
