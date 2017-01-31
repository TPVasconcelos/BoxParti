TODO:requirements.txt
Python 2.7
matplotlib
scipy
numpy

# BoxParti

*BoxParti* is a Tkinter GUI running a **_hard sphere molecular dynamics model simulation_**.

I have noticed that this version **performes better on OS X**. With minor changes it can also run on Linux and Windows. 

### If running on Windows: 
You’ll have to make some changes to the code in [BoxParti.py](BoxParti.py) to be able to open the files in the *HelpMenu*.

Around line 1483 change **_open_** to **_start_**. It should look like this: 
```python
@staticmethod
def openfile():
     # use 'start' instead of 'open' in a Windows OS
     os.system("start Python_Project_BoxParti.pdf")
```

### If running on Linux:
I have only tested this on Ubuntum distro. But the WelcomePage() class in [BoxParti.py](BoxParti.py) should only be using the **pack** or **grid** methods. This would actually be the correct way of doing it, according to Tkinter documentation!


## 4. Gallery
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


**Tomás Pereira de Vasconcelos**

tomasvasconcelos1@gmail.com 
