MultiMuC
======================
MultiMuC is a **somatic mutation caller for multi-regional tumors**.
This method can use the stochastic modeling of the state-of-the-art mutation calling methods for single-regional tumors, e.g., Strelka2, NeuSomatic, through the probability-based outputs.
Furthermore, this method can avoid a case of performance degradation ("No-TP" case defined in our paper).

How to build
----------
This script is implemented in Julia language. We checked the result in Julia v1.0.1.
Install julia language before using our script.

### Instantiate package ###
Before running the script, please instantiate the package and install dependencies.
```sh
% cd ./src
% julia
julia> ]
(v1.0) pkg>
(v1.0) pkg> activate .
(src) pkg> instantiate
```

How to run
----------
Note. At the first execution, package builging will be conducted, and simultaneous run (running in multiple processes) should be conducted after the first time of running.

To execute, type as follows.
```sh
julia --project ./src/main.jl MAP \\
--seed ${seed} -n ${iter} -b ${burnin} \\
--pdata=0.02 --pevidence=0.2 --pcon=0.999 --theta=0.5 --rho=0.1 \\
-o ${output_dir} ${BFL} ${BFH}
```
![1](/images_for_readme/readme_image_1.png)

![2](/images_for_readme/readme_image_2.png)

![3](/images_for_readme/readme_image_3.png)

![4](/images_for_readme/readme_image_4.png)

Publication
----------
The result will appear at the proceedings of ISMCO2019 (http://ismco.net/) in Springer LNCS format.

License
----------
Copyright &copy; 2019 Takuya Moriyama
Released under the [GNU General Public License, Version 3.0][GPL].
