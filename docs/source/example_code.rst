Examples
==================

Here are some examples of how to use the PyMOODENG package. These examples demonstrate the basic functionality of the package and how to set up and run simulations. All codes can be found at:

.. raw:: html

   <a href="https://colab.research.google.com/drive/1iXGNVJ6Y_ZiRfuAODnWPdClr1Cmc5zVY?usp=sharing">Colab Notebook</a>


Custom planetary system creation
--------------------------------------------------

.. code-block:: python

    from pymoodeng import plotting as p
    from pymoodeng import objects as o

    #Currently available systems
    print(p.list_systems())


    #To create a planetary system, you first define a center, around witch planets orbit
    Center = o.Body(name = "HLX-1",
                    color = "white",
                    mass = 8e20, 
                    mean_diameter = 800000)

    #Classic planet creation example
    b612 = o.Body(name = "b612",          # display name
                color = "yellow",       # display color
                mass = 1e8,             # Mass, future use
                mean_diameter = 10000,  # Size of the planet, realtive to the other bodies within the system
                orbit = o.Orbit(peri_r = 2e8, tp = 20 * 24 * 3600, e = 0.3)) #parameters: peri_r = pericentricy, tp = periodtime , e = eccentricity

    Moondeng = o.Body(name = "Moondeng",
                    color = "palevioletred", 
                    mass = 1000000, 
                    mean_diameter = 2000, 
                    orbit = o.Orbit(peri_r = 1e5, tp = 2 * 24 * 3600, e = 0.1))

    #Creation of a planetary system, which can be simulated
    #Planet - Moon 
    #name in case its further reused
    example_planetsys = o.System(name = "planetsys", # name of the system, used to reference it
            center = b612 ,   # center body 
            orbiting = Moondeng) # body or array of bodies, every object orbiting the center body

    #"Sun" - (planet - moon) 
    example_system = o.System(name = "system", 
            center = Center,
            orbiting = [example_planetsys, o.Mars] 
            )

    #Currently available systems
    print(p.list_systems())




System simulation
--------------------------------------------------

.. code-block:: python

    #DOES NOT WORK IN GOOGLE COLLAB, as it requires a pop up window to be opened
    from pymoodeng import plotting as p


    p.plot(example_system)  #Invite the created System to be animated (arguments: name (of system), END_DAY (ending day of simulation time), switching (changing between system during simulation enabled/disabled))
