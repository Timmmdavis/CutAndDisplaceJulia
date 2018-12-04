# CutAndDisplaceInJulia

#[![Build Status](https://travis-ci.com/Timmmdavis/JuliaTravisTest.svg?branch=master)](https://travis-ci.com/Timmmdavis/JuliaTravisTest)
#[![Coverage Status](https://codecov.io/gh/Timmmdavis/JuliaTravisTest/branch/master/graph/badge.svg)](https://codecov.io/gh/Timmmdavis/JuliaTravisTest)


https://github.com/KriaSoft/Folder-Structure-Conventions

/doc contains some text files on how to use this program.
/test contains the test scripts.
/src contains the Module and Function definitions.

Module dependencies are added after 'dep' in the Project.toml file:  
To add these the UID which is needed is found in the Project/Manifest.toml (directory displayed in the REPL) using `] add YYY` where YYY is the dependency. 
