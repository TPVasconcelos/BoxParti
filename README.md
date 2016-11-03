# BoxParti

**_Disclaimers:_**

**_- This is, most likely, not be up-to-date_**

**_- Runs better on OS X_** (see Section-3 if you want to run it on other systems)


## 1. What is it?
*BoxParti* is a computer program that allows users to visualise a simulation based on an ideal gas model. 

## 2. How does it work? 
The program uses the pre-installed **Python** module **Tkinter** as a *Graphical User Interface (GUI)*. This makes it very easy to display excellent plots and animations while keeping the program user-friendly and interactive.

## 3. Where can I run it?
This Python script was written to be run on **OS X**. But with minor changes it can also run on Windows and most Unix machines that have **Python 2.7** installed. Further installation of **matplotlib**, **scipy** and **numpy** modules are also required.

### If running on Linux: 
Modify WelcomePage() class on BoxParti.py and use only **pack** or only **grid** methods. This would actually be the correct way of doing it according to Tkinter documentation! If you fix this feel free to share it ;-)

### If running on Windows: 
You’ll have to make some changes to the code in [BoxParti.py](BoxParti.py) to be able to open the files in the *HelpMenu*.

In line 1483 change **_open_** to **_start_**. It should look like this: 
```python
@staticmethod
def openfile():
     # use 'start' instead of 'open' in a Windows OS
     os.system("start Python_Project_BoxParti.pdf")
```

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


## 5. FYI
This was a final project for a Scientific Computing Skills class of an Undergraduate Physics Degree. This was done when I only had been using Python for 3 months and no other real programming training. It is only the fruit of 2 to 3 weeks of work. I'm not currently developing this, but feel free to improve it! :-)


## 6. Any questions? 
Read the [Report](Python_Project_BoxParti.pdf) first!



If you have any other questions, contact me via e-mail:

**Tomás Pereira de Vasconcelos**

tomasvasconcelos1@gmail.com 
