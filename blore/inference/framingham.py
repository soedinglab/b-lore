import numpy as np

def risk(n, age, sex):
    risk = np.zeros(n)
    if age is not None and sex is not None:
        for i in range(n):
            if sex[i] == 1:
                if age[i] > 0:
                    #power = np.exp(3.06117 * (np.log(age[i]) - 3.8560))
                    #risk[i] = 1 - np.power(0.88936, power)
                    #risk[i] = 1 - np.power(0.96936, power)
                    risk[i] = 0.03 * age[i] + 0.3
            if sex[i] == 2:
                if age[i] > 0:
                    #power = np.exp(2.32888 * (np.log(age[i]) - 3.8686))
                    #risk[i] = 1 - np.power(0.95012, power)
                    #risk[i] = 1 - np.power(0.98012, power)
                    risk[i] = 0.03 * age[i]
    return risk
