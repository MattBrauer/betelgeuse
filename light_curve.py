"""The raw data comes from this page:

https://www.aavso.org/data-download

courtesy of the Royal Astronomical Society of New Zealand. I've cached it on my Github
to stop from hammering their servers."""

import requests
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

SOURCE = 'https://gist.githubusercontent.com/andyljones/112233a324d4908696750ad7a1d724e2/raw/92c91a61b52ecedba2d26962d32adc0652f8ff53/aavso_betelgeuse.txt'

# Default size is really small
plt.rcParams['figure.figsize'] = (18, 12)

def raw_data():
    path = Path('data/aavso_betelgeuse.txt')
    if not path.exists():
        print('Fetching data from Github. This is only current as of 23/12/2019!')
        path.parent.mkdir(exist_ok=True, parents=True)
        r = requests.get(SOURCE)
        r.raise_for_status()
        path.write_bytes(r.content)
    return (pd.read_csv(str(path))
                .rename(columns=lambda c: c.lower().replace(' ', '_'))
                .assign(date=lambda df: df.jd.pipe(pd.to_datetime, unit='D', origin='julian'))
                .assign(magnitude=lambda df: df.magnitude.pipe(pd.to_numeric, errors='coerce')))

def magnitudes():
    return (raw_data()
                .groupby(['date', 'observer_code']).magnitude.mean()
                .unstack())

def victors_plot(mags, start_date = '1965-01-01', end_date = '2024-01-01', window = 10):
    """Follows Victor C's methodology here: https://twitter.com/chmn_victor/status/1208898012174573568 """
    rolling = (mags
                .resample('D').mean()
                .rolling(window, min_periods=1, center=True).mean())

    well_observed = rolling[rolling.notnull().sum(1).gt(10)].reindex_like(rolling)
    inliers = well_observed.sub(well_observed.mean(1), axis=0).mean().abs().lt(.3)
    reliable = well_observed.loc[:, inliers]

    with plt.style.context('seaborn-poster'):
      ax = reliable.mean(1).plot(linewidth=1)
      ax.set_xlim(start_date, end_date)
      ax.invert_yaxis()
      ax.set_title('Betelgeuse Magnitude, Average Inlier Observation')
    #  ax.annotate('Methodology: @chmn_victor, this notebook: @andy_l_jones', (10, 10), xycoords='axes points')

if __name__ == '__main__':
    mags = magnitudes()
    victors_plot(mags)

