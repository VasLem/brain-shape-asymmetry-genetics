import pandas as pd
import numpy as np


def main():
    num_lm = 50
    num_ind = 100
    num_reps = 3
    sample_space = 100
    individual_max_deviation = 10
    side_max_deviation = 5
    rep_max_deviation = 2

    # random data points for single sample of 50 landmarks
    xyz = np.random.randint(sample_space, size=num_lm * 3)
    # generate samples from xyz, producing individual asymmetry
    noise = np.random.randint(individual_max_deviation, size=[num_ind, num_lm * 3])
    ind_xyz = pd.DataFrame(noise + xyz[np.newaxis, :])
    ind_xyz.set_index("id_" + ind_xyz.index.astype(str), inplace=True)
    # generate matching symmetry from ind_xyz, producing side asymmetry
    noise = np.random.randint(side_max_deviation, size=[1, num_lm * 3])
    match_xyz = ind_xyz + noise
    # generate fluctuating asymmetry
    noisef1 = np.random.randint(1, size=[num_ind, num_lm * 3])
    noisef2 = np.random.randint(1, size=[num_ind, num_lm * 3])
    final_side1 = ind_xyz + noisef1
    final_side1.set_index(final_side1.index.astype(str) + "_1", inplace=True)
    final_side2 = match_xyz + noisef2
    final_side2.set_index(final_side2.index.astype(str) + "_2", inplace=True)
    rep_df = pd.concat([final_side1, final_side2], axis=0)
    ret_df = None
    for cnt in range(num_reps):
        noise = np.random.randint(rep_max_deviation, size=[2 * num_ind, num_lm * 3])
        tmp = (rep_df + noise).set_index(rep_df.index.astype(str) + f"_{cnt}")
        if ret_df is None:
            ret_df = tmp
        ret_df = pd.concat([ret_df, tmp], axis=0)
    ret_df.index.name = "index"
    ret_df.columns = num_lm * ["x", "y", "z"]
    ret_df.to_csv("../SAMPLE_DATA/autogenerated.csv", sep=",")


if __name__ == "__main__":
    main()
