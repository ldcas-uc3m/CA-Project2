# Script to test the performance of the SOA program

# Note that the command "perf" should be installed: $ sudo apt install linux-tools-common
# also, change line 50 to: if(argc < 6 || (argc == 6 && argcv[6] == ">")){

echo "Running performance test for the soa program"
echo "random_seed: 2336"
echo "size_enclosure: 100000"
echo "time_step: 0.5"

# compile
g++ -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG ./sim-soa.cpp -o sim-soa 2>&1  # stderr to stdout

# run tests
for ((OBJ = 1000; OBJ<8000; OBJ = OBJ * 2)) do
    for ((ITR = 50; ITR<400; ITR = ITR * 2)) do
        for ((I = 0; I<10; I++)) do
            perf stat ./sim-soa OBJ ITR 2336 100000 0.5 >> performances.txt
            echo "\n" >> performances.txt
        done
    done
done
echo "Done"