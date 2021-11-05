# Script to test the performance of the AOS program

# Note that the command "perf" should be installed: $ sudo apt install linux-tools-common
# also, change line 49 to: if(argc < 6 || (argc == 6 && argcv[6] == ">")){

echo "Running performance test for the aos program"
echo "random_seed: 2336"
echo "size_enclosure: 100000"
echo "time_step: 0.5"

# compile
g++ -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG ./sim-aos.cpp -o sim-aos 2>&1  # stderr to stdout

# run tests
for ((OBJ = 1000; OBJ<8000; OBJ = OBJ * 2)) do
    for ((ITR = 50; ITR<400; ITR = ITR * 2)) do
        for ((I = 0; I<10; I++)) do
            sudo perf stat --output perfmance.txt ./sim-aos OBJ ITR 2336 100000 0.5

        done
    done
done
echo "Done"