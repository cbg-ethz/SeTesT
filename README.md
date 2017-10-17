# SeTesT
<p align="center">
	<img src="img/Logo.png?raw=true" alt="SeTesT"/>
</p>
David Seifert (david.seifert@bsse.ethz.ch)


## Introduction
SeTesT (**Se**lection **Tes**t in **T**ransmission) performs a test for selection in viral transmission events. It detects selection on the basis of statistically significant deviations from neutrality by accounting for strong population bottlenecks that for instance occur in the sexual transmission of HIV-1 between hosts. Standard tests such as the chi-square test and Fisher's exact test do not account for such strong population bottlenecks and yield extremely high false positive rates in practice. In addition, the aforementioned tests cannot deal with data mroe than 2 categories. SeTesT overcomes all of these limitations and can deal with various different data types.


## Idea
SeTesT is based on a 1-step Wright-Fisher process and represents a latent model. SeTesT accounts for the sampling variance inherent in next-generation sequencing. Furthermore, SeTesT also captures uncertainty due to divergent evolution occurring in the transmitter and recipient between the time of the transmission event and sampling of the viral populations.


## Features
SeTesT has the following features:

- Not prone to high false positive rates, as SeTesT takes the strong population bottleneck into account.
- Can use general phenotypes, DNA and peptide sequences as traits on which to perform a statistical test for selection.
- In the case of DNA or peptide sequences, the effects of divergent evolution can be captured using a substitution model. See the following figure as an example:
<p align="center">
	<img src="img/DivergentEvolution.png?raw=true"/>
</p>
  Here light shaded bases indicate loci that are conserved between transmitter and recipient. Darker shaded bases indicate polymorphous sites where heterogeneity occurs. Finally, gray bases indicate loci at which transmitter and recipient have diverged to a point where extinction and/or fixation of new mutant bases as occurred between the transmission event and the sampling event. SeTesT can deal with such divergent evolution to determine the most likely originating sequences in the transmitter.


## Requirements
As **SeTesT** is still under heavy development, we will not be making release tarballs yet. If you still wish to give **SeTesT** a try, you'll need the following tools and libraries:

1.  A **C++11** compliant compiler. The **SeTesT** codebase makes extensive use of C++11 features.

    GCC 4.9 and later have been verified to work, although we recommend you use at least GCC 5. Clang 3.7 and later have been verified and are also recommended, due to Clang introducing OpenMP with 3.7. Versions of Clang before 3.7 will not be able to utilise multi-threading.

2.  **Boost**; latest 1.59 release (http://www.boost.org/)

    Boost provides the necessary abstractions for many different types.

3.  **standard Unix utilities**; such as sed, etc...

    If you cannot execute a command, chances are that you are missing one of the more common utilities we require in addition to the tools listed above.

*and* either

4.  **Autoconf**; latest 2.69 release (http://www.gnu.org/software/autoconf/)

    GNU Autoconf produces the ./configure script from configure.ac.

5.  **Automake**; latest 1.15 release (http://www.gnu.org/software/automake/)

    GNU Automake produces the Makefile.in precursor, that is processed with ./configure to yield the final Makefile.

6.  **Autoconf Archive**; latest 2016.03.20 release (http://www.gnu.org/software/autoconf-archive/)

    Our configure.ac requires a number of m4 macros from the Autoconf archive.

or

4.  **CMake**; at least 3.1 (http://cmake.org)

    CMake is an alternative build system that does not require bootstrapping like the Autotools.


### OS X
We strongly recommend you use MacPorts (http://www.macports.org) to install dependencies. We also recommend you employ Clang from MacPorts, as it is the only OpenMP-capable compiler that is simultaneously ABI-compatible with installed libraries, such as boost. While building with GCC on OS X is possible, it requires an orthogonal toolchain which is far more involved and beyond the scope of this README.

### GNU/Linux
On a GNU/Linux system, the aforementioned recommendations are reversed. Most GNU/Linux distributions are built using GCC/libstdc++, which as of GCC 5.1 is not backwards compatible with Clang, and as such building with Clang produced object files will fail in the final linking step.


## Building using Autotools
1.  If you have all dependencies satisfied, proceed by generating the build system files
    ```
    ./autogen.sh
    ```

2.  Then, run the configure script. On GNU/Linux, you would do
    ```
    ./configure --prefix=<myprefix>
    ```
    where `<myprefix>` is the final root where you would like to install to.
    On OS X, you would also need to specify the OpenMP-capable C++ compiler
    ```
    ./configure --prefix=<myprefix> CXX=clang++-mp-3.7
    ```
    for instance, if you installed Clang 3.7 from MacPorts.

3.  Then, compile the sources using
    ```
    make -j2
    ```

4.  You should now have a binary called `setest` in the current build directory. You can either install this manually or call
    ```
    make install
    ```

## Building using CMake
1.  If you have all dependencies satisfied, proceed by generating creating a build dir
    ```
    mkdir build && cd build
    ```

2.  Then, run CMake
    ```
    cmake -DCMAKE_INSTALL_PREFIX=<myprefix> ..
    ```
    where `<myprefix>` is the final root where you would like to install to.
    On OS X, you would also need to specify the OpenMP-capable C++ compiler
    ```
    CXX=clang++-mp-3.7 cmake -DCMAKE_INSTALL_PREFIX=<myprefix> ..
    ```
    for instance, if you installed Clang 3.7 from MacPorts.

3.  Then, compile the sources using
    ```
    make -j2
    ```

4.  You should now have a binary called `setest` in the current build directory. You can either install this manually or call
    ```
    make install
    ```


## Running
The parameters of **SeTesT** can be viewed with the help option `-h`.
