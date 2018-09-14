# dynamicBridgeLoading
Code for the numerical calculation of dynamic bridge loading, as part of my BSc Thesis.



example_deflection_hyperloop.py: animation for bridge deflection when a hyperloop passes over

example_deflection_locomotive.py: animation for bridge deflection when a train passes over

deflection_dynamic_loading.py: fill in your own bridge and vehicle parameters at the top of 
  the script. Note that because of the type of numerical integration used (trapezoidal), 
  the model can become unstable for certain values of parameters (especially large damping 
  and spring coefficients)
  
  
  
To change the animation speed, change "interval=..." in "ani = animation.FuncAnimation(...)"
  at bottom of the script.
