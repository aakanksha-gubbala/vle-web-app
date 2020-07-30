import pandas as pd


def get_density(s, T):
    try:
        s = s.replace('%20', '+')
    except:
        pass
    density_url = 'http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=' + s
    if s == 'Hexane':
        rho = float(655)
    else:
        density = pd.read_html(density_url)[6]
        density = density.drop(density.index[0:3]).drop('No.', axis=1)
        A = float(density['A'])
        B = float(density['B'])
        C = float(density['C'])
        D = float(density['D'])
        Tmin, Tmax = float(density['Tmin']), float(density['Tmax'])  # in K

        def rho(T):
            return A / B ** (1 + (1 - T / C) ** D)
    return rho(T) if s != 'Hexane' else rho


