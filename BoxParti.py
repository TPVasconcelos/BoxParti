import os
import random
import tkinter as tk
import tkinter.messagebox as box
import warnings
from itertools import product, combinations

import matplotlib

matplotlib.use('TKAgg')

import matplotlib.animation as animation
import matplotlib.patches as patches
import numpy as np
import scipy.stats as stats
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.path import Path
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d as interp
from scipy.spatial.distance import pdist, squareform
from scipy.special import erf

warnings.filterwarnings('error')

# This fonts might not work on a Windows machine...
# you could easily find them on the internet!
TITLE_FONT = ("Helvetica", 25, "bold")
BUTTONS_FONT = ("Consolas", 25, "bold")


class BoxParti(tk.Tk):
    """
    Here we build a container where diferent frames/windows can be stacked on.
    This is not the best method because it keeps unecessary frames opened in the backgroud.
    But its good enough for this application.

    The initial frame is given by the WelcomePage class.
    New frames are raised using show_frame() func
    """

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # the container is where the frames will be stacked on
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (WelcomePage, Page2D, Page3D, HelpMenu):
            frame = F(container, self)
            self.frames[F] = frame
            # the top frame will be the visible frame
            frame.grid(row=0, column=0, sticky="nsew")

        # initial Frame will be given by the WelcomePage class
        self.show_frame(WelcomePage)

    def show_frame(self, page):
        """raises a new "page" associated to a class"""
        frame = self.frames[page]
        frame.tkraise()


class WelcomePage(tk.Frame):
    """
    The __init__ func loads the BackGround image.

    There are three 'fake' buttons on these page.
    They work by detecting the coordinates of the mouse.
    If the coordinates match the coordinates of one of the buttons, the associated frame is raised.
    I didn't find an elegant way of overcoming this issue.

    There are also two 'real' permanent buttons (QUIT & MAIN MENU)
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # Welcome BackGround Image
        welcome_image = tk.PhotoImage(file="images/BackGround.gif")
        welcome_label = tk.Label(self, image=welcome_image)
        welcome_label.image = welcome_image  # keep a reference!
        welcome_label.grid(column=0, row=0, columnspan=30, rowspan=20)

        # Detects the coordinates of the click to open the different windows
        def callback(event):
            if (1007 < event.x < 1107) and (331 < event.y < 375):
                controller.show_frame(Page2D)
            if (1007 < event.x < 1107) and (427 < event.y < 472):
                controller.show_frame(Page3D)
            if (993 < event.x < 1118) and (528 < event.y < 572):
                controller.show_frame(HelpMenu)

        welcome_label.bind("<Button-1>", callback)

        # making the two 'permanent' buttons
        button_main_menu = tk.Button(controller, text=">>>Go back to the Main Menu", font=BUTTONS_FONT,
                                     command=lambda: controller.show_frame(WelcomePage))
        button_quit = tk.Button(controller, text=">>>Quit!", command=lambda: app.destroy(),
                                font=BUTTONS_FONT)
        button_main_menu.pack(side='left')
        button_quit.pack(side='right')

        """
        # These are 'real' tk.Buttons
        button2D = tk.Button(self, text=">>> 2D Simulator",
                            command=lambda: controller.show_frame(Page_2D), font=BUTTONS_FONT)
        button3D = tk.Button(self, text=">>> 3D Simulator",
                            command=lambda: controller.show_frame(Page_3D), font=BUTTONS_FONT)
        buttonHelp = tk.Button(self, text=">>> Help Menu",
                            command=lambda: controller.show_frame(HelpMenu), font=BUTTONS_FONT)
        button2D.grid(column=29, row=10)
        button3D.grid(column=29, row=13)
        buttonHelp.grid(column=29, row=17)
        """


class Page2D(tk.Frame):
    """
    This opens the 2-Dimentional Simulator.

    Please read the discription of the individual functions for more info
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        label = tk.Label(self, text="2D Simulator", font=TITLE_FONT)
        label.pack(side="top", fill="x")

        # assign some initial conditions/values
        self.bound = [-2, 2]
        self.g = 9.81
        self.kb = 1.38e-23
        self.amu = 1.66e-27

        # assign initial lists
        self.velocity_list = []
        self.time_passed = 0
        self.trace_x = []
        self.trace_y = []

        frame = tk.Frame(self)
        # this puts all the frames as close together as possible, a bit like the tight_layout command
        frame.pack()
        # calls functions which build the 2 main program areas (Inputs and the Plot)
        self.make_plot(frame)
        self.make_inputs(frame)

    def mass(self):
        """Returns the selected mass"""
        p = self.particle.get()
        if p == "Air":
            m = 29 * self.amu
        elif p == "Oxygen":
            m = 32 * self.amu
        elif p == "Helium":
            m = 4 * self.amu
        elif p == "Tungsten Hexafluoride":
            m = 298 * self.amu
        else:  # Simulate 'air molecules' in case something goes wrong
            m = 29 * self.amu
        return m

    def make_plot(self, container):
        """making a canvas widget for the frame (placed on the left side)"""
        self.fig = Figure(figsize=(6.5, 6), dpi=100)  # makes a figure object
        # The next line is the important bit, you can not just embed a figure into a tkinter front end
        # The Canvas object can deal with loads of things but here we are going to make it
        # with our figure object embeded into it
        self.canvas = FigureCanvasTkAgg(self.fig, container)
        self.canvas.show()
        # the get_tk_widget.grid places it at co-ord (0,1) of the master frame
        self.canvas.get_tk_widget().grid(row=0, column=1)
        # adds a standard plot toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, container)
        self.toolbar.update()
        self.toolbar.grid(row=7, column=1)

    @staticmethod
    def validate(validate):
        """
        Checks if the user used the Entry boxes correctely (only floats accepted)
        validate is the input
        """
        for val in validate:
            val = str(val)
            # if nothing
            if val == "":
                box.showerror(message="Whoops! Something went wrong...\n\n"
                                      "Please remember not to leave anything in blank.\n"
                                      "Please check all your values and try again!")
                break
            for v in val:
                # accept only 0123456789 digits or '.' for the decimal place
                if v in '0123456789.':
                    try:
                        float(val)
                    except ValueError:
                        box.showerror(message="Whoops! Something went wrong...\n\n"
                                              "Please remember not to leave anything in blank.\n"
                                              "Please check all your values and try again!")
                        break
                else:
                    box.showerror(message="Whoops! Something went wrong...\n\n"
                                          "Please remember not to leave anything in blank.\n"
                                          "Please check all your values and try again!")
                    break
            else:
                continue
            break

    def make_inputs(self, container):
        """
        Here we build the frame with legends, buttons, labels and entry boxes.
        """
        # builds the frame which will hold all the inputs and their labels
        input_frame = tk.Frame(container)
        # places this frame at co-ord (0,0) of the master frame
        input_frame.grid(column=0, row=0)

        # the legend
        self.colours = tk.Label(input_frame, text="\n\nLegend:\n\n"
                                                  "YELLOW particles are the ones with speeds within one standard\n"
                                                  "deviation from the mean speed and the RED and CYAN particles\n"
                                                  "are respectively above and below one standard deviation from\n"
                                                  "the mean speed.\n\n"
                                                  "red = fast\n"
                                                  "yellow = average\n"
                                                  "cyan = slow\n", font=("Helvetica", 15, 'bold'), justify='left')
        self.colours.grid(column=0, row=13, columnspan=3)

        # the labels
        self.lblmass = tk.Label(input_frame, text="Particle=")
        self.lblmass.grid(column=0, row=0)
        self.lblN = tk.Label(input_frame, text="N=")
        self.lblN.grid(column=0, row=1)
        self.lblsize = tk.Label(input_frame, text="Size=")
        self.lblsize.grid(column=0, row=2)
        self.lblTemp = tk.Label(input_frame, text="T=")
        self.lblTemp.grid(column=0, row=3)
        self.lbldt = tk.Label(input_frame, text="dt=")
        self.lbldt.grid(column=0, row=4)
        self.lblbins = tk.Label(input_frame, text="Histogram bins=")
        self.lblbins.grid(column=0, row=5)
        self.lblv_rset = tk.Label(input_frame, text="Reset Velocities=")
        self.lblv_rset.grid(column=0, row=6)
        self.lblv_gravity = tk.Label(input_frame, text="Gravity=")
        self.lblv_gravity.grid(column=0, row=7)
        self.lblv_collisions = tk.Label(input_frame, text="Collisions=")
        self.lblv_collisions.grid(column=0, row=8)
        self.lblv_dur = tk.Label(input_frame, text="Frames=")
        self.lblv_dur.grid(column=0, row=9)

        # Entry Boxes
        self.particle = tk.StringVar(input_frame)
        self.particle.set("Air")
        self.menu_particle = tk.OptionMenu(input_frame, self.particle,
                                           "Helium", "Air", "Oxygen", "Tungsten Hexafluoride")
        self.menu_particle.grid(column=1, row=0)

        self.N = tk.Entry(input_frame)
        self.N.insert(0, "1000")
        self.N.grid(column=1, row=1)
        self.size = tk.Entry(input_frame)
        self.size.insert(1, "0.03")
        self.size.grid(column=1, row=2)
        self.temp = tk.Entry(input_frame)
        self.temp.insert(2, "300")
        self.temp.grid(column=1, row=3)
        self.dt = tk.Entry(input_frame)
        self.dt.insert(3, "0.005")
        self.dt.grid(column=1, row=4)
        self.bins = tk.Entry(input_frame)
        self.bins.insert(4, "40")
        self.bins.grid(column=1, row=5)
        self.vel_reset = tk.IntVar()
        self.vel_reset.set(1)
        self.btn_vel_reset = tk.Checkbutton(input_frame, variable=self.vel_reset)
        self.btn_vel_reset.grid(column=1, row=6)
        self.gravity = tk.IntVar()
        self.btn_gravity = tk.Checkbutton(input_frame, variable=self.gravity)
        self.btn_gravity.grid(column=1, row=7)
        self.collisions = tk.IntVar()
        self.collisions.set(1)
        self.btn_collisions = tk.Checkbutton(input_frame, variable=self.collisions)
        self.btn_collisions.grid(column=1, row=8)
        self.dur = tk.Entry(input_frame)
        self.dur.insert(6, "30")
        self.dur.grid(column=1, row=9)

        # Buttons
        self.Plot_hist = tk.Button(input_frame, text="Plot Velocities' Histogram", command=self.plot_hist_2d)
        self.Plot_hist.grid(column=0, row=10)
        self.btn_press = tk.Button(input_frame, text='Pressure', command=self.pressure_2d)
        self.btn_press.grid(column=0, row=11)
        self.btn_energy = tk.Button(input_frame, text='Kinetic Energy', command=self.energy_2d)
        self.btn_energy.grid(column=0, row=12)
        self.Plot_sim = tk.Button(input_frame, text='Start Simulator', command=self.plot_2d)
        self.Plot_sim.grid(column=1, row=10)
        self.rand_walk = tk.Button(input_frame, text='Random Walk Simulator', command=self.plot_rand_walk)
        self.rand_walk.grid(column=1, row=11)
        self.canvas_del = tk.Button(input_frame, text='Remove particles', command=lambda: self.canvas.blit())
        self.canvas_del.grid(column=1, row=12)

        self.butmass = tk.Button(input_frame, text="?",
                                 command=lambda: box.showinfo("This changes the mass of the particles",
                                                              "Helium -> 4amu\nAir -> 29amu\nOxygen -> 32amu\n"
                                                              "Tungsten Hexafluoride -> 298amu\n\n1amu = 1.66e-27g"))
        self.butmass.grid(column=2, row=0)
        self.butN = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Number of particles",
                                                                                  "If you wish to have more than 1000 particles, make sure that the 'Gravity' and "
                                                                                  "'Collisions' options are disabled. The 'Reset Velocities' option should be on.\n"
                                                                                  "This allows for faster calculations and a better animation."))
        self.butN.grid(column=2, row=1)
        self.butsize = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Size of the particles",
                                                                                     "The default value is a very exaggerated value. For an Ideal Gas simulation we would "
                                                                                     "neglect interatomic collisions."))
        self.butsize.grid(column=2, row=2)
        self.butT = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Temperature", "The temperature "
                                                                                                 "of the gas in Kelvin. In this simulation you can go down to absolute zero but note that"
                                                                                                 "the program does not simulate phase changes, so try and avoid very low temperatures."
                                                                                                 "The temperature affects the initial distribution of velocities of the particles."
                                                                                                 " (Maxwell-Boltzmann Speed Distribution).\n The simulation doesn't work for very high "
                                                                                                 "temperatures. The threshold temperature varies with the molecule/mass simulated."))
        self.butT.grid(column=2, row=3)
        self.butdt = tk.Button(input_frame, text="?",
                               command=lambda: box.showinfo("Time step", "This simulation uses a "
                                                                         "Euler Method to update the position of the particles for every dt seconds. The "
                                                                         "smaller this value is, the 'smoother' the animation will be but at the same time it "
                                                                         "will make the animation slower. For higher temperatures/higher velocities it is better "
                                                                         "to decrease dt as some particles that are supposed to undergo collision can 'pass right"
                                                                         " over each other'."))
        self.butdt.grid(column=2, row=4)
        self.butbins = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Histogram Bins", "The number of "
                                                                                                       "bins of the speed distibution histogram can be set to any value you find more "
                                                                                                       "appropriate."))
        self.butbins.grid(column=2, row=5)
        self.butrvv = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Velocity Values Reset", "When "
                                                                                                             "'ON' the program updates the values of the velocities for all the particles, deleting"
                                                                                                             " the old values. This allows the user to see the particles 'change colors' as they "
                                                                                                             "collide with each other. When displaying the Speed Distribution Histogram, the "
                                                                                                             "histogram shows the up-to-date velocities. If the button is disabled, the histogram "
                                                                                                             "will show the accumulation of all velocities since the beginning of the simulation."))
        self.butrvv.grid(column=2, row=6)
        self.butgrav = tk.Button(input_frame, text="?",
                                 command=lambda: box.showinfo("Gravity", "When 'ON', the gravity "
                                                                         "is exaggerated by a factor of 100 for visualisation purposes. Not essential for the "
                                                                         "simulation of an Ideal Gas."))
        self.butgrav.grid(column=2, row=7)
        self.butcoli = tk.Button(input_frame, text="?",
                                 command=lambda: box.showinfo("Collisions", "When 'ON', particles"
                                                                            " can undergo collisions with each other. This should be always 'ON' when visualizing"
                                                                            " the 'Random Walk' option."))
        self.butcoli.grid(column=2, row=8)
        self.butdur = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Frames", "Number of frames in the "
                                                                                              "animation/simulation. This is an unfortunate feature that had to be implemented in "
                                                                                              "this program. This forces the animation to crash after n frames."))
        self.butdur.grid(column=2, row=9)

    def mb_speed_dist(self, v):
        """ returns Maxwell-Boltzmann speed distribution for speeds """
        T = float(self.temp.get())
        kb = self.kb
        m = self.mass()
        return (m / (2 * np.pi * kb * T)) ** 1.5 * 4 * np.pi * v ** 2 * np.exp(-m * v ** 2 / (2 * kb * T))

    def mb_cdf(self):
        """returns CMF of the Maxwell-Boltzmann speed distribution"""
        # validade temperature entry
        self.validate([self.temp.get()])

        T = float(self.temp.get())
        v = np.arange(0, 50000, 1)
        kb = self.kb
        m = self.mass()
        a = np.sqrt(kb * T / m)
        return erf(v / (np.sqrt(2) * a)) - np.sqrt(2 / np.pi) * v * np.exp(-v ** 2 / (2 * a ** 2)) / a

    def vel_generator(self):
        """Returns and array of velocities in 2D [vx, vy]"""
        # create Inverse of the CDF first by interpolating it
        v = np.arange(0, 50000, 1)
        cdf = self.mb_cdf()
        inv_cdf = interp(cdf, v)

        # Generator
        N = int(self.N.get())
        rand_nums = np.random.random((N, 2))
        speeds = inv_cdf(rand_nums)

        new_speeds = np.empty((2,))
        for n in range(0, N):
            """ This gives equally distributed random speeds in the positive and negative directions"""
            vx = (speeds[n, 0] / 100.) * random.choice([-1, 1])  # scaling factor of 10^-2
            vy = (speeds[n, 1] / 100.) * random.choice([-1, 1])  # scaling factor of 10^-2
            vel_vec = np.array([vx, vy])
            new_speeds = np.vstack((new_speeds, vel_vec))
        new_speeds = np.delete(new_speeds, 0, axis=0)

        # the mean value should be close to zero since the velocities
        # are equaly distrybited into positive and negative values
        # mean_val = mean(speeds)
        # std_val = np.std(speeds)

        return new_speeds

    def update_2d(self, dt):
        """
        2-Dimensional Update

        update_2d(dt)
        type(dt) = float or int

        Updates the position of N particles taking in consideration boundary
        conditions and perfect elastic collisions.

        Returns an array with the positions of N particles
        And stores the updated velocities of N particles
        """
        # validate entries
        self.validate([self.N.get(), self.size.get(), self.dur.get()])
        # get size and number of particels
        size = float(self.size.get())
        N = int(self.N.get())

        if self.time_passed == 0:
            # Set up initial state/conditions (random coords and MB speeds)
            coord = np.random.random((N, 2)) * 4. - 2.
            speeds = self.vel_generator()
            self.state = np.hstack((coord, speeds))

        while self.time_passed < float(self.dur.get()):
            self.time_passed += 1  # update time
            # ignore numpy Runtime and Warning errors
            with np.errstate(invalid='ignore', divide='ignore', under='ignore'):

                ### update position ###
                if self.gravity.get():
                    self.state[:, :2] += dt * self.state[:, 2:]
                    # Gravity exagerated by a factor of 1000 for visualisation purposes
                    self.state[:, 3] += - 0.5 * (self.g * 1000) * (dt ** 2)
                else:
                    self.state[:, :2] += dt * self.state[:, 2:]

                ### reflect at boundary ###
                # if TRUE: put the particles back at the boundary
                # AND invert the velocity component
                # at x = -b
                test = (self.state[:, 0] - size < self.bound[0])
                self.state[test, 0] = self.bound[0] + size
                self.state[test, 2] *= -1
                # at x = b
                test = (self.state[:, 0] + size > self.bound[1])
                self.state[test, 0] = self.bound[1] - size
                self.state[test, 2] *= -1
                # at y = -b
                test = (self.state[:, 1] - size < self.bound[0])
                self.state[test, 1] = self.bound[0] + size
                self.state[test, 3] *= -1
                # at y = b
                test = (self.state[:, 1] + size > self.bound[1])
                self.state[test, 1] = self.bound[1] - size
                self.state[test, 3] *= -1

                ### elastic collisions  ###
                if self.collisions.get():  # if colisions 'activated'
                    # reletive distance between all points
                    D = squareform(pdist(self.state[:, :2]))
                    # D = np.nan_to_num(D)
                    # check if they are undergoing a colision
                    d_where1, d_where2 = np.where(D < 2 * size)
                    unique = (d_where1 < d_where2)
                    d_where1 = d_where1[unique]
                    d_where2 = d_where2[unique]
                    coli = zip(d_where1, d_where2)
                    # update velocities of colliding pairs
                    for i, j in coli:
                        # mass
                        m1 = self.mass()
                        m2 = self.mass()
                        # position vectors
                        r1 = self.state[i, :2]
                        r2 = self.state[j, :2]
                        # velocity vectors
                        v1 = self.state[i, 2:]
                        v2 = self.state[j, 2:]
                        # relative position and velocity vectors
                        r_rel = r1 - r2
                        v_rel = v1 - v2
                        # momentum vector of the center of mass
                        p1 = m1 * v1
                        p2 = m2 * v2
                        v_cm = (p1 + p2) / (m1 + m2)
                        # collisions of perfect elastic hard spheres
                        rr_rel = np.dot(r_rel, r_rel)
                        vr_rel = np.dot(v_rel, r_rel)
                        v_rel = 2 * r_rel * vr_rel / rr_rel - v_rel
                        # assign new velocity vectors
                        self.state[i, 2:] = v_cm + v_rel * m2 / (m1 + m2)
                        self.state[j, 2:] = v_cm - v_rel * m1 / (m1 + m2)

            ### get net velocities ###
            if self.vel_reset.get():
                # ressets velocity list
                self.velocity_list = []
            for n in range(0, N):
                # append new speeds
                new_vel = np.sqrt(self.state[n, 2] ** 2 + self.state[n, 3] ** 2)
                self.velocity_list.append(new_vel)

            ### and Finaly... ###
            return self.state[:, :2]

        self.time_passed = 0  # reset clock to zero
        self.trace_x = []  # needed for the random walk simulator
        self.trace_y = []  # needed for the random walk simulator
        self.fig.clf()
        # we need to force the animation to crash/finish
        # Still need to find an elegant solution to this issue
        # ...ideally, a 'STOP' button
        raise

    @staticmethod
    def init_2d():
        """animation init function"""
        return []

    def animate_2d(self, i):
        """Animation Functon for the 2D simulator"""
        self.validate([self.dt.get()])

        dt = float(self.dt.get())
        data = self.update_2d(dt)[:, :2]

        # create empty lists
        data_very_fast = []
        data_very_fast_x = []
        data_very_fast_y = []
        data_normal = []
        data_normal_x = []
        data_normal_y = []
        data_very_slow = []
        data_very_slow_x = []
        data_very_slow_y = []

        # get statistical values
        mean_val = np.mean(self.velocity_list)
        std_val = np.std(self.velocity_list)
        fast = mean_val + 1.5 * std_val
        slow = mean_val - 1.5 * std_val

        # separate the different particles into different velocities
        size = int(self.N.get())
        for v in range(0, size):
            if self.velocity_list[v] > fast:
                data_very_fast.append([data[v, 0], data[v, 1]])
            elif self.velocity_list[v] > slow:
                data_normal.append([data[v, 0], data[v, 1]])
            else:
                data_very_slow.append([data[v, 0], data[v, 1]])

        # separate list into two (for x and y)
        for n in range(0, len(data_very_fast)):
            data_very_fast_x.append(data_very_fast[n][0])
            data_very_fast_y.append(data_very_fast[n][1])
        for n in range(0, len(data_normal)):
            data_normal_x.append(data_normal[n][0])
            data_normal_y.append(data_normal[n][1])
        for n in range(0, len(data_very_slow)):
            data_very_slow_x.append(data_very_slow[n][0])
            data_very_slow_y.append(data_very_slow[n][1])

        return ax.plot(data_very_fast_x, data_very_fast_y, 'ro',
                       data_normal_x, data_normal_y, 'yo',
                       data_very_slow_x, data_very_slow_y, 'co',
                       ms=3, markeredgewidth=0.0)

    def plot_2d(self):
        """Plot 2D animation"""
        global ax
        ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
                                  xlim=(-0.2 + self.bound[0], 0.2 + self.bound[1]),
                                  ylim=(-0.2 + self.bound[0], 0.2 + self.bound[1]))

        # Draw a frame for the boundaries
        verts = [(self.bound[0], self.bound[0]),  # left, bottom
                 (self.bound[0], self.bound[1]),  # left, top
                 (self.bound[1], self.bound[1]),  # right, top
                 (self.bound[1], self.bound[0]),  # right, bottom
                 (self.bound[0], self.bound[0]),  # ignored
                 ]
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, fc='none', lw=1)
        ax.add_patch(patch)

        anim = animation.FuncAnimation(self.fig, self.animate_2d, interval=1,
                                       blit=True, init_func=self.init_2d)
        self.canvas.show()

    def animate_hist_2d(self, i):
        """Animation Function for the speed distribution Histogram"""
        self.validate([self.temp.get(), self.dt.get(), self.bins.get()])

        dt = float(self.dt.get())
        self.update_2d(dt)
        T = float(self.temp.get())
        # transform nan floats to zeros
        data = np.nan_to_num(self.velocity_list)

        ## Using MB_Speed_Dist(v) func to plot the expected distribution
        # v = np.arange(0, max(data)*100, 1)
        # fv = self.MB_Speed_Dist(v)

        # scipy's stats func maxwell to fit data
        maxwell = stats.maxwell
        x = np.linspace(0, max(data), 10000)
        params = maxwell.fit(data)

        ax = self.fig.add_subplot(111)
        ax.clear()
        ax.plot(x, maxwell.pdf(x, *params), "r-", label='T=' + str(T) + 'K', lw=2)
        ax.hist(data, bins=int(self.bins.get()), normed=1)

        ax.set_xlim(left=0, right=25)
        ax.set_ylim(bottom=0, top=0.55)
        ax.tick_params(labelsize='small')
        ax.set_xlabel('$v (m/s)$')
        ax.set_ylabel('$f(v)/(m/s)$')
        ax.legend(loc=0)

    def plot_hist_2d(self):
        """Plot speed distribution Histogram"""
        anim = animation.FuncAnimation(self.fig, self.animate_hist_2d, interval=1)
        self.canvas.show()  # updates the canvas

    def animate_rand_walk(self, i):
        """random walk animation function"""
        self.validate([self.dt.get()])

        dt = float(self.dt.get())
        data = self.update_2d(dt)[:, :2]

        # Draw trace
        self.trace_x.append(data[0, 0])
        self.trace_y.append(data[0, 1])
        t = ax.plot(self.trace_x, self.trace_y, "b-", lw=1.3, alpha=0.5)

        # Plot new position
        d = ax.plot(data[0, 0], data[0, 1], 'bo', ms=9, markeredgewidth=0.0)

        return d + t

    def plot_rand_walk(self):
        """plot random walk"""
        global ax

        ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
                                  xlim=(-0.2 + self.bound[0], 0.2 + self.bound[1]),
                                  ylim=(-0.2 + self.bound[0], 0.2 + self.bound[1]))

        # Draw a frame for the boundaries
        verts = [(self.bound[0], self.bound[0]),  # left, bottom
                 (self.bound[0], self.bound[1]),  # left, top
                 (self.bound[1], self.bound[1]),  # right, top
                 (self.bound[1], self.bound[0]),  # right, bottom
                 (self.bound[0], self.bound[0]),  # ignored
                 ]
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, fc='none', lw=1)
        ax.add_patch(patch)

        anim = animation.FuncAnimation(self.fig, self.animate_rand_walk, interval=1,
                                       blit=True, init_func=self.init_2d)
        self.canvas.show()

    def pressure_2d(self):
        """pop-up pressure measurements"""
        try:
            self.validate([self.N.get(), self.temp.get()])

            N = int(self.N.get())
            N = float(N)
            T = int(self.temp.get())
            V = abs((2 * self.bound[1]) ** 3)  # Volume
            m = self.mass()
            kb = self.kb
            R = 8.314

            # turn list into a numpy array
            new_vel_list = np.array(self.velocity_list)
            # mean squared velocity
            vm = (1. / 2) * np.mean((new_vel_list * 100) ** 2)

            p_vel = (1. / 3) * (N / V) * m * vm
            p_temp = (N * kb * T) / V
            error = abs(p_temp - p_vel)

            box.showinfo("Pressure", "Computed from the average velocity\n(Kinetic Theory):\n-->{0}Pa"
                                     "\n\nComputed from the Temperature\n(Ideal Gas Law):\n-->{1}Pa"
                                     "\n\nERROR: {2}".format(p_vel, p_temp, error))

        except RuntimeWarning:
            box.showwarning("Warning!", "Please run the simulation at least once before getting reading of "
                                        "Pressure or Kinetic Energy.")

    def energy_2d(self):
        """pop-up kinetic energy measurments"""
        try:
            self.validate([self.N.get(), self.temp.get()])

            N = int(self.N.get())
            N = float(N)
            T = int(self.temp.get())
            V = abs(200 * self.bound[1] ** 3)  # Volume
            m = self.mass()
            kb = self.kb

            # turn list into a numpy array
            new_vel_list = np.array(self.velocity_list)
            vm = (1. / 2) * np.mean((new_vel_list * 100) ** 2)

            e_vel = (1. / 2) * m * vm
            e_temp = (3. / 2) * kb * T
            e_total = e_temp * N
            error = abs(e_vel - e_temp)

            box.showinfo("Translational Kinetic Energy", "Computed from the average velocity\n(Kinetic Theory):"
                                                         "\n-->{0}J"
                                                         "\n\nComputed from the Temperature\n(Ideal Gas Law):\n-->{1}J"
                                                         "\n\nERROR: {2}"
                                                         "\n\nTotal Translational Kinetic Energy:\n-->{3}J".format(
                e_vel, e_temp, error, e_total))

        except RuntimeWarning:
            box.showwarning("Warning!", "Please run the simulation at least once before getting reading of "
                                        "Pressure or Kinetic Energy.")


class Page3D(tk.Frame):
    """
    This class opens the 3-Dimentional Simulator.

    Please read the discription of the individual functions for more info.

    NOTE: Page_2D class has more detailed comments!
          Please refer to those, as most of this code is the same but for 3 Dimentions!

    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        label = tk.Label(self, text="3D Simulator", font=TITLE_FONT)
        label.pack(side="top", fill="x")

        # assign some initial conditions/values
        self.bound = [-2, 2]
        self.g = 9.81
        self.kb = 1.38e-23
        self.amu = 1.66e-27

        # assign initial lists
        self.velocity_list = []
        self.time_passed = 0
        self.trace_x = []
        self.trace_y = []
        self.trace_z = []

        frame = tk.Frame(self)
        # this puts all the frames as close together as possible, a bit like the tight_layout command
        frame.pack()
        # calls functions which build the 2 main program areas
        self.make_plot(frame)
        self.make_inputs(frame)

    def mass(self):
        """Returns the selected mass"""
        p = self.particle.get()
        if p == "Air":
            m = 29 * self.amu
        elif p == "Oxygen":
            m = 32 * self.amu
        elif p == "Helium":
            m = 4 * self.amu
        elif p == "Tungsten Hexafluoride":
            m = 297.83 * self.amu
        else:  # Simulate 'air molecules' in case something goes wrong
            m = 29 * self.amu
        return m

    def make_plot(self, container):
        """making a canvas widget for the frame (placed on the left side)"""
        # makes a canvas widget that will show the final calculation
        self.fig = Figure(figsize=(6.5, 6), dpi=100)  # makes a figure object
        # The next line is the important bit, you can not just embed a figure into a tkinter front end
        # The Canvas object can deal with loads of things but here we are going to make it
        # with our figure object embeded into it
        self.canvas = FigureCanvasTkAgg(self.fig, container)
        self.canvas.show()
        # the get_tk_widget.grid places it at co-ord (0,1) of the master frame
        self.canvas.get_tk_widget().grid(row=0, column=1)
        # adds a standard plot toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, container)
        self.toolbar.update()
        self.toolbar.grid(row=7, column=1)

    @staticmethod
    def validate(validate):
        """Checks if the user user the Entry boxes correctely with only floats accepted"""
        for val in validate:
            val = str(val)
            if val == "":
                box.showerror(message="Whoops! Something went wrong...\n\n"
                                      "Please remember not to leave anything in blank.\n"
                                      "Please check all your values and try again!")
                break
            for v in val:
                if v in '0123456789.':
                    try:
                        float(val)
                    except ValueError:
                        box.showerror(message="Whoops! Something went wrong...\n\n"
                                              "Please remember not to leave anything in blank.\n"
                                              "Please check all your values and try again!")
                        break
                else:
                    box.showerror(message="Whoops! Something went wrong...\n\n"
                                          "Please remember not to leave anything in blank.\n"
                                          "Please check all your values and try again!")
                    break
            else:
                continue
            break

    def make_inputs(self, container):
        """
        Here we build the frame with legends, buttons, labels and entry boxes.
        """
        # builds the frame which will hold all the inputs and their labels
        input_frame = tk.Frame(container)
        # places this frame at co-ord (0,0) of the master frame
        input_frame.grid(column=0, row=0)

        self.colours = tk.Label(input_frame, text="\n\nLegend:\n\n"
                                                  "YELLOW particles are the ones with speeds within one standard\n"
                                                  "deviation from the mean speed and the RED and CYAN particles\n"
                                                  "are respectively above and below one standard deviation from\n"
                                                  "the mean speed.\n\n"
                                                  "red = fast\n"
                                                  "yellow = average\n"
                                                  "cyan = slow\n", font=("Helvetica", 15, 'bold'), justify='left')
        self.colours.grid(column=0, row=13, columnspan=3)

        # the lables
        self.lblmass = tk.Label(input_frame, text="Particle=")
        self.lblmass.grid(column=0, row=0)
        self.lblN = tk.Label(input_frame, text="N=")
        self.lblN.grid(column=0, row=1)
        self.lblsize = tk.Label(input_frame, text="Size=")
        self.lblsize.grid(column=0, row=2)
        self.lblTemp = tk.Label(input_frame, text="T=")
        self.lblTemp.grid(column=0, row=3)
        self.lbldt = tk.Label(input_frame, text="dt=")
        self.lbldt.grid(column=0, row=4)
        self.lblbins = tk.Label(input_frame, text="Histogram bins=")
        self.lblbins.grid(column=0, row=5)
        self.lblv_rset = tk.Label(input_frame, text="Reset Velocities=")
        self.lblv_rset.grid(column=0, row=6)
        self.lblv_gravity = tk.Label(input_frame, text="Gravity=")
        self.lblv_gravity.grid(column=0, row=7)
        self.lblv_collisions = tk.Label(input_frame, text="Collisions=")
        self.lblv_collisions.grid(column=0, row=8)
        self.lblv_dur = tk.Label(input_frame, text="Frames=")
        self.lblv_dur.grid(column=0, row=9)

        # Entry Boxes
        self.particle = tk.StringVar(input_frame)
        self.particle.set("Air")
        self.menu_particle = tk.OptionMenu(input_frame, self.particle,
                                           "Helium", "Air", "Oxygen", "Tungsten Hexafluoride")
        self.menu_particle.grid(column=1, row=0)
        self.N = tk.Entry(input_frame)
        self.N.insert(0, "1000")
        self.N.grid(column=1, row=1)
        self.size = tk.Entry(input_frame)
        self.size.insert(1, "0.03")
        self.size.grid(column=1, row=2)
        self.temp = tk.Entry(input_frame)
        self.temp.insert(2, "300")
        self.temp.grid(column=1, row=3)
        self.dt = tk.Entry(input_frame)
        self.dt.insert(3, "0.005")
        self.dt.grid(column=1, row=4)
        self.bins = tk.Entry(input_frame)
        self.bins.insert(4, "40")
        self.bins.grid(column=1, row=5)
        self.vel_reset = tk.IntVar()
        self.vel_reset.set(1)
        self.btn_vel_reset = tk.Checkbutton(input_frame, variable=self.vel_reset)
        self.btn_vel_reset.grid(column=1, row=6)
        self.gravity = tk.IntVar()
        self.btn_gravity = tk.Checkbutton(input_frame, variable=self.gravity)
        self.btn_gravity.grid(column=1, row=7)
        self.collisions = tk.IntVar()
        self.collisions.set(1)
        self.btn_collisions = tk.Checkbutton(input_frame, variable=self.collisions)
        self.btn_collisions.grid(column=1, row=8)
        self.dur = tk.Entry(input_frame)
        self.dur.insert(6, "30")
        self.dur.grid(column=1, row=9)

        # Buttons
        self.Plot_hist = tk.Button(input_frame, text="Plot Velocities' Histogram", command=self.plot_hist_3d)
        self.Plot_hist.grid(column=0, row=10)
        self.btn_press = tk.Button(input_frame, text='Pressure', command=self.pressure_3d)
        self.btn_press.grid(column=0, row=11)
        self.btn_energy = tk.Button(input_frame, text='Kinetic Energy', command=self.energy_3d)
        self.btn_energy.grid(column=0, row=12)
        self.Plot_sim = tk.Button(input_frame, text='Start Simulator', command=self.plot_3d)
        self.Plot_sim.grid(column=1, row=10)
        self.rand_walk = tk.Button(input_frame, text='Random Walk Simulator', command=self.plot_rand_walk)
        self.rand_walk.grid(column=1, row=11)
        self.canvas_del = tk.Button(input_frame, text='Remove particles', command=lambda: self.canvas.blit())
        self.canvas_del.grid(column=1, row=12)

        self.butmass = tk.Button(input_frame, text="?",
                                 command=lambda: box.showinfo("This changes the mass of the particles",
                                                              "Helium -> 4amu\nAir -> 29amu\nOxygen -> 32amu\n"
                                                              "Tungsten Hexafluoride -> 298amu\n\n1amu = 1.66e-27g"))
        self.butmass.grid(column=2, row=0)
        self.butN = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Number of particles",
                                                                                  "If you wish to have more than 1000 particles, make sure that the 'Gravity' and "
                                                                                  "'Collisions' options are disabled. The 'Reset Velocities' option should be on.\n"
                                                                                  "This allows for faster calculations and a better animation."))
        self.butN.grid(column=2, row=1)
        self.butsize = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Size of the particles",
                                                                                     "The default value is a very exaggerated value. For an Ideal Gas simulation we would "
                                                                                     "neglect interatomic collisions."))
        self.butsize.grid(column=2, row=2)
        self.butT = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Temperature", "The temperature "
                                                                                                 "of the gas in Kelvin. In this simulation you can go down to absolute zero but note that"
                                                                                                 "the program does not simulate phase changes, so try and avoid very low temperatures."
                                                                                                 "The temperature affects the initial distribution of velocities of the particles."
                                                                                                 " (Maxwell-Boltzmann Speed Distribution).\n The simulation doesn't work for very high "
                                                                                                 "temperatures. The threshold temperature varies with the molecule/mass simulated."))
        self.butT.grid(column=2, row=3)
        self.butdt = tk.Button(input_frame, text="?",
                               command=lambda: box.showinfo("Time step", "This simulation uses a "
                                                                         "Euler Method to update the position of the particles for every dt seconds. The "
                                                                         "smaller this value is, the 'smoother' the animation will be but at the same time it "
                                                                         "will make the animation slower. For higher temperatures/higher velocities it is better "
                                                                         "to decrease dt as some particles that are supposed to undergo collision can 'pass right"
                                                                         " over each other'."))
        self.butdt.grid(column=2, row=4)
        self.butbins = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Histogram Bins", "The number of "
                                                                                                       "bins of the speed distibution histogram can be set to any value you find more "
                                                                                                       "appropriate."))
        self.butbins.grid(column=2, row=5)
        self.butrvv = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Velocity Values Reset", "When "
                                                                                                             "'ON' the program updates the values of the velocities for all the particles, deleting"
                                                                                                             " the old values. This allows the user to see the particles 'change colors' as they "
                                                                                                             "collide with each other. When displaying the Speed Distribution Histogram, the "
                                                                                                             "histogram shows the up-to-date velocities. If the button is disabled, the histogram "
                                                                                                             "will show the accumulation of all velocities since the beginning of the simulation."))
        self.butrvv.grid(column=2, row=6)
        self.butgrav = tk.Button(input_frame, text="?",
                                 command=lambda: box.showinfo("Gravity", "When 'ON', the gravity "
                                                                         "is exaggerated by a factor of 100 for visualisation purposes. Not essential for the "
                                                                         "simulation of an Ideal Gas."))
        self.butgrav.grid(column=2, row=7)
        self.butcoli = tk.Button(input_frame, text="?",
                                 command=lambda: box.showinfo("Collisions", "When 'ON', particles"
                                                                            " can undergo collisions with each other. This should be always 'ON' when visualizing"
                                                                            " the 'Random Walk' option."))
        self.butcoli.grid(column=2, row=8)
        self.butdur = tk.Button(input_frame, text="?", command=lambda: box.showinfo("Frames", "Number of frames in the "
                                                                                              "animation/simulation. This is an unfortunate feature that had to be implemented in "
                                                                                              "this program. This forces the animation to crash after n frames."))
        self.butdur.grid(column=2, row=9)

    def mb_speed_dist(self, v):
        """ Maxwell-Boltzmann speed distribution for speeds """
        T = float(self.temp.get())
        kb = self.kb
        m = self.mass()
        return 100 * (m / (2 * np.pi * kb * T)) ** 1.5 * 4 * np.pi * v ** 2 * np.exp(-m * v ** 2 / (2 * kb * T))

    def mb_cdf(self):
        """Cumulative Distribution function of the Maxwell-Boltzmann speed distribution"""
        self.validate([self.temp.get()])

        T = float(self.temp.get())
        v = np.arange(0, 2500, 1)
        kb = self.kb
        m = self.mass()
        a = np.sqrt(kb * T / m)
        return erf(v / (np.sqrt(2) * a)) - np.sqrt(2 / np.pi) * v * np.exp(-v ** 2 / (2 * a ** 2)) / a

    def vel_generator(self):
        """Returns and array of velocities in 2D [vx, vy]"""
        # create Inverse of the CDF first by interpolating it
        v = np.arange(0, 2500, 1)
        cdf = self.mb_cdf()
        inv_cdf = interp(cdf, v)

        # Generator
        N = int(self.N.get())
        rand_nums = np.random.random((N, 3))
        speeds = inv_cdf(rand_nums)

        new_speeds = np.empty((3,))
        for n in range(0, N):
            a = (speeds[n, 0] / 100.) * random.choice([-1, 1])  # scaling factor of 10^-2
            b = (speeds[n, 1] / 100.) * random.choice([-1, 1])  # scaling factor of 10^-2
            c = (speeds[n, 2] / 100.) * random.choice([-1, 1])  # scaling factor of 10^-2
            vel_vec = np.array([a, b, c])
            new_speeds = np.vstack((new_speeds, vel_vec))
        new_speeds = np.delete(new_speeds, 0, axis=0)

        # the mean value should be close to zero since the velocities
        # are equaly distrybited into positive and negative values
        # mean_val = mean(speeds)
        # std_val = np.std(speeds)

        return new_speeds

    def update_3d(self, dt):
        """
        3-Dimensional Update

        update_3d(dt)
        type(dt) = float or int

        Updates the position of N particles taking in consideration boundary
        conditions and perfect elastic collisions.

        Returns an array with the positions of N particles
        And stores the updated velocities of N particles
        """
        self.validate([self.N.get(), self.size.get(), self.dur.get()])

        size = float(self.size.get())
        N = int(self.N.get())
        if self.time_passed == 0:
            coord = np.random.random((N, 3)) * 4. - 2.
            speeds = self.vel_generator()
            self.state = np.hstack((coord, speeds))

        while self.time_passed < float(self.dur.get()):
            self.time_passed += 1
            with np.errstate(invalid='ignore', divide='ignore', under='ignore'):
                ### update position ###
                if self.gravity.get():
                    self.state[:, :3] += dt * self.state[:, 3:]
                    # Gravity exagerated by a factor of 1000 for visual porpuses
                    self.state[:, 5] += - 0.5 * (self.g * 1000) * (dt ** 2)
                else:
                    self.state[:, :3] += dt * self.state[:, 3:]

                ### reflect at boundary ###
                # at x = -b
                test = (self.state[:, 0] - size < self.bound[0])
                self.state[test, 0] = self.bound[0] + size
                self.state[test, 3] *= -1
                # at x = b
                test = (self.state[:, 0] + size > self.bound[1])
                self.state[test, 0] = self.bound[1] - size
                self.state[test, 3] *= -1

                # at y = -b
                test = (self.state[:, 1] - size < self.bound[0])
                self.state[test, 1] = self.bound[0] + size
                self.state[test, 4] *= -1
                # at y = b
                test = (self.state[:, 1] + size > self.bound[1])
                self.state[test, 1] = self.bound[1] - size
                self.state[test, 4] *= -1

                # at z = -b
                test = (self.state[:, 2] - size < self.bound[0])
                self.state[test, 2] = self.bound[0] + size
                self.state[test, 5] *= -1
                # at z = b
                test = (self.state[:, 2] + size > self.bound[1])
                self.state[test, 2] = self.bound[1] - size
                self.state[test, 5] *= -1

                ### elastic collisions ###
                if self.collisions.get():
                    # find pairs of particles undergoing a collision
                    D = squareform(pdist(self.state[:, :3]))
                    # D = np.nan_to_num(D)
                    d_where1, d_where2 = np.where(D < 2 * size)
                    unique = (d_where1 < d_where2)
                    d_where1 = d_where1[unique]
                    d_where2 = d_where2[unique]
                    coli = zip(d_where1, d_where2)
                    # update velocities of colliding pairs
                    for i, j in coli:
                        # mass
                        m1 = self.mass()
                        m2 = self.mass()
                        # location vector
                        r1 = self.state[i, :3]
                        r2 = self.state[j, :3]
                        # velocity vector
                        v1 = self.state[i, 3:]
                        v2 = self.state[j, 3:]
                        # relative location & velocity vectors
                        r_rel = r1 - r2
                        v_rel = v1 - v2
                        # momentum vector of the center of mass
                        p1 = m1 * v1
                        p2 = m2 * v2
                        v_cm = (p1 + p2) / (m1 + m2)
                        # collisions of hard-spheres
                        # reflect v_rel over r_rel
                        rr_rel = np.dot(r_rel, r_rel)
                        vr_rel = np.dot(v_rel, r_rel)
                        v_rel = 2 * r_rel * vr_rel / rr_rel - v_rel
                        # assign new velocities
                        self.state[i, 3:] = v_cm + v_rel * m2 / (m1 + m2)
                        self.state[j, 3:] = v_cm - v_rel * m1 / (m1 + m2)

            ### get net velocities ###
            if self.vel_reset.get():
                self.velocity_list = []
            for n in range(0, N):
                new_vel = np.sqrt(self.state[n, 3] ** 2 + self.state[n, 4] ** 2 + self.state[n, 5] ** 2)
                self.velocity_list.append(new_vel)

            ### Finalize ###
            return self.state[:, :3]

        self.time_passed = 0
        self.trace_x = []
        self.trace_y = []
        self.trace_z = []
        self.fig.clf()
        raise  # force the animation to crash

    @staticmethod
    def init_3d():
        """animation init function"""
        return []

    def animate_3d(self, i):
        """Animation Functon for the 2D simulator"""
        self.validate([self.dt.get()])

        dt = float(self.dt.get())
        data = self.update_3d(dt)[:, :3]

        data_very_fast = []
        data_very_fast_x = []
        data_very_fast_y = []
        data_very_fast_z = []
        data_normal = []
        data_normal_x = []
        data_normal_y = []
        data_normal_z = []
        data_very_slow = []
        data_very_slow_x = []
        data_very_slow_y = []
        data_very_slow_z = []

        mean_val = np.mean(self.velocity_list)
        std_val = np.std(self.velocity_list)
        fast = mean_val + 1.5 * std_val
        slow = mean_val - 1.5 * std_val
        size = int(self.N.get())
        for v in range(0, size):
            if self.velocity_list[v] > fast:
                data_very_fast.append([data[v, 0], data[v, 1], data[v, 2]])
            elif self.velocity_list[v] > slow:
                data_normal.append([data[v, 0], data[v, 1], data[v, 2]])
            else:
                data_very_slow.append([data[v, 0], data[v, 1], data[v, 2]])

        for n in range(0, len(data_very_fast)):
            data_very_fast_x.append(data_very_fast[n][0])
            data_very_fast_y.append(data_very_fast[n][1])
            data_very_fast_z.append(data_very_fast[n][2])
        for n in range(0, len(data_normal)):
            data_normal_x.append(data_normal[n][0])
            data_normal_y.append(data_normal[n][1])
            data_normal_z.append(data_normal[n][2])
        for n in range(0, len(data_very_slow)):
            data_very_slow_x.append(data_very_slow[n][0])
            data_very_slow_y.append(data_very_slow[n][1])
            data_very_slow_z.append(data_very_slow[n][2])

        sf = ax.plot(data_very_fast_x, data_very_fast_y, data_very_fast_z, 'ro', ms=5, markeredgewidth=0.0)
        n = ax.plot(data_normal_x, data_normal_y, data_normal_z, 'yo', ms=5, markeredgewidth=0.0)
        ss = ax.plot(data_very_slow_x, data_very_slow_y, data_very_slow_z, 'co', ms=5, markeredgewidth=0.0)

        all_3d_pts = sf + n + ss

        return all_3d_pts

    def plot_3d(self):
        """Plot 2D animation"""
        global ax
        ax = Axes3D(self.fig)

        # draw cube
        r = [self.bound[0], self.bound[1]]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s - e)) == r[1] - r[0]:
                ax.plot3D(*zip(s, e), color="k", lw=1.5)

        ax.set_xlim3d([-0.2 + self.bound[0], 0.2 + self.bound[1]])
        ax.set_xlabel('X')
        ax.set_ylim3d([-0.2 + self.bound[0], 0.2 + self.bound[1]])
        ax.set_ylabel('Y')
        ax.set_zlim3d([-0.2 + self.bound[0], 0.2 + self.bound[1]])
        ax.set_zlabel('Z')

        def cube_faces():
            """ This is for presentation purposes only
                If you wish to NOT draw faces of a cube 'comment' Cube_Faces() below """
            # Face 1
            x1 = np.array([[-2, -2, -2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            y1 = np.array([[-2, 2, 2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            z1 = np.array([[-2, -2, 2, 2, -2],
                           [-2, -2, -2, -2, -2]])
            # Face 2
            x2 = np.array([[-2, 2, 2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            y2 = np.array([[2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2]])
            z2 = np.array([[-2, -2, 2, 2, -2],
                           [-2, -2, -2, -2, -2]])
            # Face 3
            x3 = np.array([[-2, -2, 2, 2, -2],
                           [-2, -2, -2, -2, -2]])
            y3 = np.array([[-2, 2, 2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            z3 = np.array([[-2, -2, -2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            ax.plot_surface(x1, y1, z1, color="k", alpha=0.6)
            ax.plot_surface(x2, y2, z2, color="k", alpha=0.4)
            ax.plot_surface(x3, y3, z3, color="k", alpha=0.2)

        cube_faces()

        anim = animation.FuncAnimation(self.fig, self.animate_3d, interval=1, blit=True, init_func=self.init_3d)
        self.canvas.show()

    def animate_hist_3d(self, i):
        """Animation Function for the speed distribution Histogram"""
        self.validate([self.temp.get(), self.dt.get(), self.bins.get()])

        dt = float(self.dt.get())
        self.update_3d(dt)
        T = float(self.temp.get())

        data = np.nan_to_num(self.velocity_list)

        maxwell = stats.maxwell
        x = np.linspace(0, max(data), 3000)
        params = maxwell.fit(data)

        ax = self.fig.add_subplot(111)
        ax.clear()
        ax.plot(x, maxwell.pdf(x, *params), "r-", label='T=' + str(T) + 'K', lw=2)
        ax.hist(data, bins=int(self.bins.get()), normed=1)

        ax.set_xlim(left=0, right=20)
        ax.set_ylim(bottom=0, top=0.4)
        ax.tick_params(labelsize='small')
        ax.set_xlabel('$v (m/s)$')
        ax.set_ylabel('$f(v)/(m/s)$')
        ax.legend(loc=0)

    def plot_hist_3d(self):
        """Plot speed distribution Histogram"""
        anim = animation.FuncAnimation(self.fig, self.animate_hist_3d, interval=1, repeat=False)
        self.canvas.show()  # updates the canvas,

    def animate_rand_walk(self, i):
        """random walk animation function"""
        self.validate([self.dt.get()])

        dt = float(self.dt.get())
        data = self.update_3d(dt)[:, :3]

        # Draw trace
        self.trace_x.append(data[0, 0])
        self.trace_y.append(data[0, 1])
        self.trace_z.append(data[0, 2])
        t = ax.plot(self.trace_x, self.trace_y, self.trace_z, "b-", lw=1.5, alpha=0.6)

        # Plot new position

        d = ax.plot([data[0, 0]], [data[0, 1]], [data[0, 2]], 'bo', ms=8, markeredgewidth=1)

        p = d + t
        return p

    def plot_rand_walk(self):
        """plot random walk"""
        global ax
        ax = Axes3D(self.fig)

        # draw cube
        r = [self.bound[0], self.bound[1]]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s - e)) == r[1] - r[0]:
                ax.plot3D(*zip(s, e), color="g", lw=1)

        ax.set_xlim3d([-0.2 + self.bound[0], 0.2 + self.bound[1]])
        ax.set_xlabel('X')
        ax.set_ylim3d([-0.2 + self.bound[0], 0.2 + self.bound[1]])
        ax.set_ylabel('Y')
        ax.set_zlim3d([-0.2 + self.bound[0], 0.2 + self.bound[1]])
        ax.set_zlabel('Z')

        def cube_faces():
            """ This is for presentation purposes only
                If you wish to NOT draw faces of a cube 'comment' Cube_Faces() below """
            # Face 1
            x1 = np.array([[-2, -2, -2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            y1 = np.array([[-2, 2, 2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            z1 = np.array([[-2, -2, 2, 2, -2],
                           [-2, -2, -2, -2, -2]])
            # Face 2
            x2 = np.array([[-2, 2, 2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            y2 = np.array([[2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2]])
            z2 = np.array([[-2, -2, 2, 2, -2],
                           [-2, -2, -2, -2, -2]])
            # Face 3
            x3 = np.array([[-2, -2, 2, 2, -2],
                           [-2, -2, -2, -2, -2]])
            y3 = np.array([[-2, 2, 2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            z3 = np.array([[-2, -2, -2, -2, -2],
                           [-2, -2, -2, -2, -2]])
            ax.plot_surface(x1, y1, z1, color="k", alpha=0.6)
            ax.plot_surface(x2, y2, z2, color="k", alpha=0.4)
            ax.plot_surface(x3, y3, z3, color="k", alpha=0.2)

        cube_faces()

        anim = animation.FuncAnimation(self.fig, self.animate_rand_walk, interval=1,
                                       blit=True, init_func=self.init_3d)
        self.canvas.show()

    def pressure_3d(self):
        """pop-up pressure measurements"""
        try:
            self.validate([self.N.get(), self.temp.get()])

            N = int(self.N.get())
            N = float(N)
            T = int(self.temp.get())
            V = abs((2 * self.bound[1]) ** 3)
            m = self.mass()
            kb = self.kb
            R = 8.314

            new_vel_list = np.array(self.velocity_list)
            vm = (1. / 3) * np.mean((new_vel_list * 100) ** 2)

            p_vel = (1. / 3) * (N / V) * m * vm
            p_temp = (N * kb * T) / V
            error = abs(p_vel - p_temp)

            box.showinfo("Pressure", "Computed from the average velocity\n(Kinetic Theory):\n-->{0}Pa"
                                     "\n\nComputed from the Temperature\n(Ideal Gas Law):\n-->{1}Pa"
                                     "\n\nERROR: {2}".format(p_vel, p_temp, error))

        except RuntimeWarning:
            box.showwarning("Warning!", "Please run the simulation at least once before getting reading of "
                                        "Pressure or Kinetic Energy.")

    def energy_3d(self):
        """pop-up kinetic energy measurments"""
        try:
            self.validate([self.N.get(), self.temp.get()])

            N = int(self.N.get())
            N = float(N)
            T = int(self.temp.get())
            V = abs(200 * self.bound[1] ** 3)
            m = self.mass()
            kb = self.kb

            new_vel_list = np.array(self.velocity_list)
            vm = (1. / 3) * np.mean((new_vel_list * 100) ** 2)

            e_vel = (1. / 2) * m * vm
            e_temp = (3. / 2) * kb * T
            e_total = e_temp * N
            error = abs(e_vel - e_temp)

            box.showinfo("Translational Kinetic Energy", "Computed from the average velocity\n(Kinetic Theory):"
                                                         "\n-->{0}J"
                                                         "\n\nComputed from the Temperature\n(Ideal Gas Law):\n-->{1}J"
                                                         "\n\nERROR: {2}"
                                                         "\n\nTotal Translational Kinetic Energy:\n-->{3}J".format(
                e_vel, e_temp, error, e_total))

        except RuntimeWarning:
            box.showwarning("Warning!", "Please run the simulation at least once before getting reading of "
                                        "Pressure or Kinetic Energy.")


class HelpMenu(tk.Frame):
    """
    The HelpMenu Frame class is not very developed.
    From here, you have access to a YouTube video showing how the program works and a small writen Report.
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # pack the UnderDevelopment.gif image
        dev_image = tk.PhotoImage(file="images/UnderDevelopment.gif")
        dev_label = tk.Label(self, image=dev_image)
        dev_label.image = dev_image  # keep a reference!
        dev_label.pack()

        label = tk.Label(self, text="Sorry! This Page is still under development...\n\nIf you're having trouble with "
                                    "any of the simulations,\nrefer to the 'help buttons' located next to each option."
                                    "\n\nIf the problem persists open this document:", font=TITLE_FONT)
        label.pack(side="top", fill="x")

        frame = tk.Frame(self)
        # this puts all the frames as close together as possible, a bit like the tight_layout command
        frame.pack()
        # calls functions which build the 2 main program areas
        self.make_inputs(frame)

    def make_inputs(self, container):
        input_frame = tk.Frame(container)
        input_frame.pack()

        # Report Button
        report_btn = tk.Button(input_frame, text=">>> This button directs you to the Project Report", font=BUTTONS_FONT,
                               command=self.openfile)
        report_btn.pack(side="top", pady=10)

        # Video Example Button
        video_image = tk.PhotoImage(file="images/movie.gif")
        video_btn = tk.Button(input_frame, compound="right", text="This is a Video Example (YouTube link) -->",
                              image=video_image, command=self.openvideo)
        video_btn.image = video_image  # it's important keep a 'reference'!!
        video_btn.pack(side="top", pady=40)

    @staticmethod
    def openfile():
        import platform
        if platform.system().lower() in ('windows', 'win32'):
            os.system("start Python_Project_BoxParti.pdf")
        else:
            os.system("open Python_Project_BoxParti.pdf")

    @staticmethod
    def openvideo():
        import webbrowser
        webbrowser.open("https://youtu.be/iL4_eIdm02E")
        # If you have the video in your computer and dont have access to the internet:
        # use 'start' instead of 'open' in a Windows OS
        # os.system("open ExampleVideo.mov")


if __name__ == "__main__":
    app = BoxParti()
    app.mainloop()
