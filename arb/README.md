# Arb numerics

These are some initial work by Brad from 2020/2021, where he was extending the
'special functions' approach of Devreese, in order to try and evaluate
frequency dependent mobility functions, at the level of the counter integration
in FHIP. 

Shortly after all this work, he realised that actually brute-forcing the
intengral with QuadGK *BEFORE* making the shift to a contour integration had
far better numerics, and was faster and had simpler code than using these
special functions (because we had to take them to ridiculous accuracy). 

These special functions supplement the Julia ArbNumerics package, but we never
pushed the code upstream. :(

