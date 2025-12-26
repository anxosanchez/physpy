# data/database.py
SOLVENTS = {
    "Agua": {
        "MW": 18.015, "Tc": 647.1, "Pc": 220.6, "Omega": 0.344, "Z_RA": 0.2338,
        "rho_ref": 997.0, "visc_ref": 0.89, "sigma_ref": 72.8, "Paracor": 52.0,
        "dD": 15.5, "dP": 16.0, "dH": 42.3
    },
    "Etanol": {
        "MW": 46.069, "Tc": 513.9, "Pc": 61.4, "Omega": 0.644, "Z_RA": 0.2520,
        "rho_ref": 789.0, "visc_ref": 1.07, "sigma_ref": 22.1, "Paracor": 128.0,
        "dD": 15.8, "dP": 8.8, "dH": 19.4
    },
    "Acetona": {
        "MW": 58.08, "Tc": 508.1, "Pc": 47.0, "Omega": 0.304, "Z_RA": 0.2547,
        "rho_ref": 784.0, "visc_ref": 0.31, "sigma_ref": 23.3, "Paracor": 161.0,
        "dD": 15.5, "dP": 10.4, "dH": 7.0
    },
    "Tolueno": {
        "MW": 92.14, "Tc": 591.8, "Pc": 41.0, "Omega": 0.264, "Z_RA": 0.2646,
        "rho_ref": 867.0, "visc_ref": 0.56, "sigma_ref": 28.5, "Paracor": 246.0,
        "dD": 18.0, "dP": 1.4, "dH": 2.0
    },
    "Metanol": {
        "MW": 32.04, "Tc": 512.6, "Pc": 80.9, "Omega": 0.556, "Z_RA": 0.2329,
        "rho_ref": 791.0, "visc_ref": 0.54, "sigma_ref": 22.5, "Paracor": 91.0,
        "dD": 15.1, "dP": 12.3, "dH": 22.3
    },
    "n-Butanol": {
        "MW": 74.12, "Tc": 563.0, "Pc": 44.2, "Omega": 0.594, "Z_RA": 0.2587,
        "rho_ref": 810.0, "visc_ref": 2.54, "sigma_ref": 24.6, "Paracor": 210.0,
        "dD": 16.0, "dP": 5.7, "dH": 15.8
    },
    "Acetato de Butilo": {
        "MW": 116.16, "Tc": 579.0, "Pc": 31.1, "Omega": 0.434, "Z_RA": 0.2590,
        "rho_ref": 881.0, "visc_ref": 0.73, "sigma_ref": 25.2, "Paracor": 315.0,
        "dD": 15.8, "dP": 3.7, "dH": 6.3
    },
    "MEK": {
        "MW": 72.11, "Tc": 535.5, "Pc": 41.5, "Omega": 0.323, "Z_RA": 0.2600,
        "rho_ref": 805.0, "visc_ref": 0.41, "sigma_ref": 24.0, "Paracor": 198.0,
        "dD": 16.0, "dP": 9.0, "dH": 5.1
    },
    "Xileno": {
        "MW": 106.16, "Tc": 617.0, "Pc": 35.1, "Omega": 0.302, "Z_RA": 0.2630,
        "rho_ref": 861.0, "visc_ref": 0.62, "sigma_ref": 28.7, "Paracor": 284.0,
        "dD": 17.6, "dP": 1.0, "dH": 3.1
    },
    "Isopropanol": {
        "MW": 60.10, "Tc": 508.3, "Pc": 47.6, "Omega": 0.665, "Z_RA": 0.2501,
        "rho_ref": 785.0, "visc_ref": 2.04, "sigma_ref": 21.7, "Paracor": 160.0,
        "dD": 15.8, "dP": 6.1, "dH": 16.4
    },
    "Ciclohexano": {
        "MW": 84.16, "Tc": 553.5, "Pc": 40.7, "Omega": 0.210, "Z_RA": 0.2729,
        "rho_ref": 778.1, "visc_ref": 0.89, "sigma_ref": 24.3, "Paracor": 217.0,
        "dD": 16.8, "dP": 0.0, "dH": 0.2
    },
    "n-Heptano": {
        "MW": 100.20, "Tc": 540.2, "Pc": 27.4, "Omega": 0.349, "Z_RA": 0.2635,
        "rho_ref": 684.0, "visc_ref": 0.39, "sigma_ref": 19.7, "Paracor": 311.0,
        "dD": 15.3, "dP": 0.0, "dH": 0.0
    },
    "Acetato de Etilo": {
        "MW": 88.11, "Tc": 523.3, "Pc": 38.8, "Omega": 0.366, "Z_RA": 0.2554,
        "rho_ref": 900.3, "visc_ref": 0.42, "sigma_ref": 23.4, "Paracor": 216.0,
        "dD": 15.8, "dP": 5.3, "dH": 7.2
    },
    "THF": {
        "MW": 72.11, "Tc": 540.2, "Pc": 51.9, "Omega": 0.226, "Z_RA": 0.2560,
        "rho_ref": 889.0, "visc_ref": 0.46, "sigma_ref": 26.4, "Paracor": 174.0,
        "dD": 16.8, "dP": 5.7, "dH": 8.0
    },
    "DMSO": {
        "MW": 78.13, "Tc": 719.0, "Pc": 56.5, "Omega": 0.345, "Z_RA": 0.2350,
        "rho_ref": 1100.0, "visc_ref": 1.99, "sigma_ref": 42.9, "Paracor": 189.0,
        "dD": 18.4, "dP": 16.4, "dH": 10.2
    }
}

