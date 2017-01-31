# DFSPH

Divergence-Free SPH for simulation of water based on the method presented by J. Bender & D. Koschier: [https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf](https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf). Made in the course TN1008 year 2016.

The project was built in visual studio 2015 in c++ using the built in compiler. 

To run the simulation a few external librarys need to be included:

* Since the simulation use opengl, glew and glfw these are needed.
* The GLM library is used for some caluclations. 
* [CompactNSearch](https://github.com/InteractiveComputerGraphics/CompactNSearch) is used for neighboor search.
* [Imgui](https://github.com/ocornut/imgui) was used to create the interface.


Video:

[![DFSPH youtube video](http://img.youtube.com/vi/7Vduk3ByZug/0.jpg)](http://www.youtube.com/watch?v=7Vduk3ByZug)
