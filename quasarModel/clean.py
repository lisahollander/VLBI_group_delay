import numpy as np

def analytical_visibility(df, U, V, header):
    c_inv = 1.0 / (3 * 10**8)
    dTheta = (np.multiply.outer(U, df['DELTAX'].values*np.pi/180) + np.multiply.outer(V, df['DELTAY'].values*np.pi/180)) * c_inv
    flux = df['FLUX'].values

    theta = dTheta * header.w
    cs, sn = flux * np.cos(theta), flux * np.sin(theta)
    return np.sum(cs, axis=2), -np.sum(sn, axis=2), -np.sum(sn * dTheta, axis=2), -np.sum(cs * dTheta, axis=2)

def get_clean_results(df,header, U,V):
    grid = np.zeros((header.size, header.size))  # Use np.nan or 0 as an initial value

    # Populate the grid with values from 'z' at corresponding (x, y) positions
    for _ , row in df.iterrows():
        grid[int(row['Y']), int(row['X'])] = row['FLUX']

    reV, imV, dReV, dImV = analytical_visibility(df,U, V,header)

    clean_visibility = reV + 1j*imV
    dclean_visibility = 1/(reV**2 + imV**2)*(dImV*reV - dReV*imV)

    return np.fliplr(clean_visibility), np.fliplr(dclean_visibility)

