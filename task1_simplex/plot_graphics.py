import matplotlib.pyplot as plt


def main():

    num_of_threads = [1, 2, 4, 8]

    simplex_basic_omp_50000 = [1.00578006, 1.833656, 3.2292131, 4.01454]
    simplex_basic_omp_100000 = [0.996102, 1.91879092, 3.28196929, 4.35677224]
    simplex_basic_omp_1000000 = [0.99028, 1.874465, 3.209318, 4.06627]

    plt.plot(num_of_threads, simplex_basic_omp_50000)
    plt.plot(num_of_threads, simplex_basic_omp_100000)
    plt.plot(num_of_threads, simplex_basic_omp_1000000)
    plt.legend(["50000 iterations", "100000 iterations", "1000000 iterations"], loc="best")
    plt.savefig("simplex_basic.png", dpi = 90)
    plt.show()

    simplex_modified_omp_50000 = [1.842015, 4.338867, 5.4644261, 6.781684]
    simplex_modified_omp_100000 = [1.832739, 3.351518, 5.35788864, 6.85931996]
    simplex_modified_omp_1000000 = [1.85726, 3.393428, 5.530534, 7.06112969]

    plt.plot(num_of_threads, simplex_modified_omp_50000)
    plt.plot(num_of_threads, simplex_modified_omp_100000)
    plt.plot(num_of_threads, simplex_modified_omp_1000000)
    plt.legend(["50000 iterations", "100000 iterations", "1000000 iterations"], loc="best")
    plt.savefig("simplex_modified.png", dpi = 90)
    plt.show()


if __name__ == "__main__":
    main()