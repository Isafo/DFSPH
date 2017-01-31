# DFSPH

Divergence-Free SPH for simulation of water based on the method presented in: [https://animation.rwth-aachen.de/media/papers/2016-TVCG-ViscousDFSPH.pdf](https://animation.rwth-aachen.de/media/papers/2016-TVCG-ViscousDFSPH.pdf). Made in the course TN1008 year 2016.

The project was built in visual studio 15 in c++ using the built in compiler. 

To run the simulation a few external librarys need to be included.

1 Since the simulation use opengl, glew and glfw these are needed.
2 The GLM library is used for some caluclations. 
3 CompactNSearch is used for neighboor search.
4 Imgui was used to create a interface
