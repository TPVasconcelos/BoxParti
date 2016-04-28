# BoxParti
*This is open source!*

*Please do use this and modify it as you wish.*

*If you do use it and modify it, I would appreciate if you shared your work with me!*

*And let me know if you can find a way of overcoming some of the issues in the code:*

- *Once the animation starts running on tkinter, I could not find a way of stoping it with a ’STOP’ button. I do not like the fact that I have to ask the user to chose the animation duration before hand (using the ‘frames’ option).*
- *For some reason there is a backend conflict when I ran this on Canopy in OS X. However, I can run it on the command line or on PyCharm with no problem (I do need to specify the 'TKAgg’ backend)*
- *In the 3D Simulation I would like to make it possible to zoom in and out and move around the plot.*
- *I know it’s bad practice using global function, as it increases the chance of getting a bug on the program. But I couldn’t get around not using it on functions:*
	- *Page2D.plot_2d 		(global ax)*
	- *Page2D.plot_rand_walk 	(global ax)*
	- *Page3D.plot_3d 		(global ax)*
	- *Page3D.plot_rand_walk 	(global ax)*

## 1. What is it?
*BoxParti* is a computer program that allows users to visualise a simulation based on an ideal gas model. 

## 2. How does it work? 
The program uses the pre-installed **Python** module **Tkinter** as a *Graphical User Interface (GUI)*. This makes it very easy to display excellent plots and animations while keeping the program user-friendly and interactive.

## 3. Where can I run it?
This Python script was written to be run on **OS X**. But with minor changes it can also run on Windows and most Unix machines that have **Python 2.7** installed. Further installation of the most up-to-date versions of **matplotlib**, **scipy** and **numpy** modules may also be required for the script to run with no errors.
### If running on windows: 
You’ll have to make some changes to the code in [BoxParti.py](BoxParti.py) to be able to open the files in the *HelpMenu*.

In line 1483 change **_open_** to **_start_**. It should look like this: 
```python
@staticmethod
def openfile():
     # use 'start' instead of 'open' in a Windows OS
     os.system("start Python_Project_BoxParti.pdf")
```
**_NOTE: this might not be up-to-date!_**

## 4. Gallery
### The HomePage/WelcomePage
![welcome page](readme_gallery/WelcomePage.jpg)

### The 2D Simulator Page
![2D Page](readme_gallery/2DPage.jpg)

### 2D Simulation with collisions
![2D](readme_gallery/2D_SIM.gif)

### 3D Simulation with collisions and gravity
![3D](readme_gallery/3D_SIM.gif)

### 2D Brownian Motion
![Brownian Motion](readme_gallery/brownian.gif)

### Maxwell-Boltzmann Speed Distribution
![MB](readme_gallery/MB.gif)


## 5. Preamble… 
I had only been using Python for 3 months when I did this project.

I worked on it for around 2 to 3 weeks.

This was a final project for a Scientific Computing Skills class of a  Second Year Physics Degree.


## 6. Any questions? 
Read the [Report](Python_Project_BoxParti.pdf) first!



If you have any other questions, contact me via e-mail:

Tomás Pereira de Vasconcelos

tomasvasconcelos1@gmail.com 
