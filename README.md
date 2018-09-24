**Deprecated** - I have expanded this project into [mdsea](https://github.com/TPVasconcelos/mdsea), a stand-alone Python molecular dynamics library equipped with multiple analysis and visualisation tools.

---

# BoxParti

_"My first ever programming project."_ - A Tkinter GUI running a hard-sphere molecular dynamics simulation.


I developed this project for a "Scientific Computing Skills" module in the first year of my Physics Undergraduate degree. I have not modified it since the module's deadline. I have only corrected some minor bugs, to make it easier for others to try it out and play with the source.

This was developed on macOS and I've noticed that the GUI does not render so well on Linux and Windows machines. This could be fixed and will be _"left as an exercise for the reader"_. \[**Hint:** The `WelcomePage()` class in [BoxParti.py](BoxParti.py) should only use the **pack** or **grid** methods (not both). This would actually be the correct way of doing things, [according to the Tkinter documentation](http://effbot.org/tkinterbook/grid.htm).\]

If you're a beginner to Python programming, feel free to fork this project and play around with the source, to see some examples of tkinter, numpy, and scipy libraries and Python's OOP. Unless you find a **minor** bug or improvement which does not affect the overall project and source, please do not make any pull requests, as this project is no longer being developed.


## Some snapshots of the interface and simulations:

#### The HomePage/WelcomePage
![welcome page](readme_gallery/WelcomePage.jpg)

#### The 2D Simulator Page
![2D Page](readme_gallery/2DPage.jpg)

#### 2D Simulation with collisions
![2D](readme_gallery/2D_SIM.gif)

#### 3D Simulation with collisions and gravity
![3D](readme_gallery/3D_SIM.gif)

#### 2D Brownian Motion
![Brownian Motion](readme_gallery/brownian.gif)

#### Maxwell-Boltzmann Speed Distribution
![MB](readme_gallery/MB.gif)
