import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np 
import vtolParam as P


class vtolAnimation:
    '''
        Create vtol animation
    '''
    def __init__(self):
        self.flagInit = True                  # Used to indicate initialization
        self.fig, self.ax = plt.subplots()    # Initializes a figure and axes object
        self.handle = []                      # Initializes a list object that will
                                              # be used to contain handles to the
                                              # patches and line objects.
        self.cl = P.cl                        # center pod length
        self.ch = P.ch                        # center pod altitude
        self.d = P.d                          # distance between center pods and rotors
        self.h0 = P.h0                        # initial vertical position of center pod
        plt.axis([-10*self.cl, 10*self.cl, -0.1, 10*self.h0])  # Change the x,y axis limits
        plt.plot([-10*self.cl, 10*self.cl], [0, 0], 'b--')     # Draw a base line
        plt.xlabel('z')
        plt.ylabel('h')
        # plt.axis('equal')

        # Draw vtol is the main function that will call the functions:
        # drawCenterPod, drawRightRotor, drawLeftRotor, and drawConnectingRods to create the animation.
    def drawVtol(self, u):
        # Process inputs to function
        z = u[0]        # Horizontal position of pod, m
        h = u[1]        # Vertical position of pod
        theta = u[2]    # Angle of pendulum, rads

        self.drawCenterPod(z, h, theta)
        self.drawRightRotor(z, h, theta)
        self.drawLeftRotor(z, h, theta)
        self.drawConnectingRods(z, h, theta)
        self.ax.axis('equal') # This will cause the image to not distort

        # After each function has been called, initialization is over.
        if self.flagInit == True:
            self.flagInit = False

    def drawCenterPod(self, z, h, theta):
        X = [z - self.cl/2*np.cos(theta), z + self.cl/2*np.cos(theta)]   # X data points
        Y = [h + self.cl/2*np.sin(theta), h - self.cl/2*np.sin(theta)]  # Y data points

        # When the class is initialized, a line object will be
        # created and added to the axes. After initialization, the
        # line object will only be updated.
        if self.flagInit == True:
            # Create the line object and append its handle
            # to the handle list.
            line, = self.ax.plot(X, Y, lw=15, c='blue')
            self.handle.append(line)
        else:
            self.handle[0].set_xdata(X)  # Update the line
            self.handle[0].set_ydata(Y)

    def drawRightRotor(self, z, h, theta):
        x = z + (self.cl/2 + self.d)*np.cos(theta)   # x coordinate
        y = h - self.cl/2*np.sin(theta) - self.d*np.sin(theta)   # y coordinate
        xy = (x,y)                                   # Center of circle

        # When the class is initialized, a CirclePolygon patch object will
        # be created and added to the axes. After initialization, the
        # CirclePolygon patch object will only be updated.
        if self.flagInit == True:
            # Create the CirclePolygon patch and append its handle
            # to the handle list
            self.handle.append(mpatches.Ellipse(xy,
              width=0.2, height=0.1,
              angle=-theta * 180/np.pi,
              fc='limegreen', ec='black'))
            self.ax.add_patch(self.handle[1])  # Add the patch to the axes
        else:
            self.handle[1].center = xy
            self.handle[1].angle = -theta * 180/np.pi

    def drawLeftRotor(self, z, h, theta):
        x = z - (self.cl/2 + self.d)*np.cos(theta)   # x coordinate
        y = h + self.cl/2*np.sin(theta) + self.d*np.sin(theta)   # y coordinate
        xy = (x,y)                                   # Center of circle

        # When the class is initialized, a CirclePolygon patch object will
        # be created and added to the axes. After initialization, the
        # CirclePolygon patch object will only be updated.
        if self.flagInit == True:
            # Create the CirclePolygon patch and append its handle
            # to the handle list
            self.handle.append(mpatches.Ellipse(xy,
              width=0.2, height=0.1,
              angle=-theta * 180/np.pi,
              fc='limegreen', ec='black'))
            self.ax.add_patch(self.handle[2])  # Add the patch to the axes
        else:
            self.handle[2].center = xy
            self.handle[2].angle = -theta * 180/np.pi

    def drawConnectingRods(self, z, h, theta):
        XR = [z + self.cl/2*np.cos(theta), z + (self.cl/2 + self.d)*np.cos(theta)]   # X data points
        YR = [h - self.cl/2*np.sin(theta), h - self.cl/2*np.sin(theta) - self.d*np.sin(theta)]  # Y data points
        XL = [z - self.cl / 2 * np.cos(theta), z - (self.cl / 2 + self.d) * np.cos(theta)]  # X data points
        YL = [h + self.cl / 2 * np.sin(theta), h + self.cl / 2 * np.sin(theta) + self.d*np.sin(theta)]  # Y data points

        # When the class is initialized, a line object will be
        # created and added to the axes. After initialization, the
        # line object will only be updated.
        if self.flagInit == True:
            # Create the line object and append its handle
            # to the handle list.
            line_r, = self.ax.plot(XR, YR, lw=2, c='black')
            self.handle.append(line_r)
            line_l, = self.ax.plot(XL, YL, lw=2, c='black')
            self.handle.append(line_l)
        else:
            self.handle[3].set_xdata(XR)  # Update the line
            self.handle[3].set_ydata(YR)
            self.handle[4].set_xdata(XL)  # Update the line
            self.handle[4].set_ydata(YL)

# Used see the animation from the command line
if __name__ == "__main__":

    simAnimation = vtolAnimation()    # Create Animate object
    z = 0.0                               # Position of cart, m
    theta = 0.0*np.pi/180               # Angle of pendulum, rads
    h = 0                               # Vtol altitude
    d = 1                               # Distance between center pod and rotors
    simAnimation.drawVtol([z, h, theta])  # Draw the Vtol
    #plt.show()
    # Keeps the program from closing until the user presses a button.
    print('Press key to close')
    plt.waitforbuttonpress()
    plt.close()