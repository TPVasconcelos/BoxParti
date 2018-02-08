**Deprecated** - I have expanded this project into [pymodys](https://github.com/TPVasconcelos/pymodys), a python molecular dynamics simulation engine that also allows for 2D (matplotlib) and 3D (Blender) animations.

---

# BoxParti

_My first programming project._

This was a project for a "Scientific Computing Skills" module in the first year of my BSc Physics. It is a Tkinter GUI running a **hard sphere molecular dynamics model simulation**.

I have noticed that this version performes better on macOS. With minor changes it can also run smoothly on Linux and Windows machines. 


**If running on Linux** - I have only tested this on Ubuntum. But the `WelcomePage()` class in [BoxParti.py](BoxParti.py) should be using only either the **pack** or **grid** methods (not both). This would actually be the correct way of doing things, according to the Tkinter documentation.


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
