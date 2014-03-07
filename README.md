This readme was created by Jason to try to give some description of the code for the reader, and for himself!

However, you should be careful. I may forget to update this if I make a change!

Note: ZrCuAl2011.eam.alloy is the correct potential file.

Since our method is statistical, we really need multiple HRMC simulations so that we can compare them. We want to get our simulations as "different" as possible, and hope that the final structures are similar. One method of doing this is starting with a crystalline and an amorphous structure as the initial model.


##Rotations
The number of rotations in the model can vary, we typically set it to 1 40 20 = 211 rotations. That means the intensity calculation is being done that many times, once for each rotation, and then combined into the V(k) curve for experimental comparison.

Note that the first rotation should be 0 0 0, so that the first model in the rotation array (mrot) is identical to the model in rmc.f90. Yes, thats a tad bit of a memory waste but really no big deal.

##Pixels
Ideally, we want every atom in the model to contribute once and only once to the intensity (per rotation). For RMC this is *essential*. With round pixels we will either have overlapping atoms (square inscribed in circle idea) or we will leave out atoms (cirle inscribed in square). By the way, the model is a box (it cant be a sphere due to boundary effects) and thats where the squares come from.

So for RMC we use square pixels to remove this effect, even though it does not exactly replicate the experiment. Removing this effect for RMC is essential because if an atom moves into a space that is not in the intensity calculation, or is duplicated in the intensity calculation, the intensity is per atom can be seen to be different. This is not good. Note the best explanation...

For Femsim, we do want to exactly replicate the experiment and we are simply calculating the variance of the model. Therefore, we are not moving atoms (and have no input data by the way) and so the RMC probelms are irrelevant here. We will use overlapping pixels instead of leaving atoms out here.

In the code we want the pixel placement to be as obvious for the user as possible. Page 11 in my notebook shows pictures of the pixel setup for a 3x3 case.

There is something I dont like about the pixel size / placement in the model currently. Right now (as well as before my modifications), for a square pixel, the pixels are placed like they are bigger than they actually are. We set the size to res*sqrt(2), but this is the length of the diagonal of a square, not the actual side length of a square. Whether we have a round pixel or a square pixel, the size should be the same because the radius of the circle and the side length of the square are the same. Instead, we are writing things as if the square side length were bigger than it really is. Its ok, because the model size was increased to res*sqrt(2) and the intensity calculation includes all atoms within the res*sqrt(2) larger square, but this really doesnt make sense. We should set the size of the square pixel the same as res because it is the same. However, then all our model sizes would need to be changed to the smaller square (which means less atoms). Everything works now, it just isnt super logical unfortunately. I could be wrong here - this might have something to do with why we cant pick our model size for femsim, I dont understand why we cant right now.

Now, Paul actually wants to have large overlap between the pixels so that every part of the model gets sampled multiple times (4?) (for femsim - NOT rmc). Note though that the atoms should be sampled ~ the same number of times! And we use round pixels in femsim, not square.

However, for RMC we pick the model size and we want each atom to be sampled once and only once. Because we pick the model size, we make it so that the model size is an integer multiple of the resolution (pixel size). This way the pixels fit in perfectly. Right now the code does this (with the exception of the above comment).
If you allow pixels to go outside the model then you need to change how things are done elsewhere as well due to PBCs which are not enabled for the hutch_list_pixels functions.

##Alpha
Alpha is the fitting parameter or weighting titled "!weight of fem; alpha" in the param_file.in.

The higher alpha is, the more the V(k) will be fitted / taken into account in the cost function (compared to the energy of the model).

We have been trying to figure out a good number to set it at. I was using 550 for some of the simulations, but I still think this is a fairly significant overfit to the V(k). We believe a better number is 0.0256eV/atom*#atoms. 0.0256eV/atom was picked because it is kT at room temp, which is about how much energy an atom has at RT compared to 0 Kelvin (right?). For a model with, for example, 1523 atoms, alpha=39.0. Accuracy is fairly insignificant here most likely.


## Acceptance rate
We ideally want an acceptance rate around 50% throughout the HRMC simulation. This is especially difficult to acheive at low temps. We can reduce the move size for lower temps, which we do to some extent already. The acceptance rate also starts near 60-70%.

I can consider making the max move size a function of the acceptance rate. That way I can try to keep the acceptance rate at 50%. However, super small moves dont do a helluva lot of good either sometimes. They can though.

##Energy
The final energy of the model after an HRMC simulation, which may take 3 million steps (attempted atom moves), is up for debate. Whats important is the energy of the model compared to the energy of a model of the same composition, but minimized against only the energy (LAMMPS may be usable for this). We hope that the energy will converge within 0.0256eV/atom, which is what alpha is based on. Even better is better. Perfect agreement is best. JWHs models were about 0.018eV/atom higher than the perfect if Paul remembers correctly.

The energy calcuation is not the same as LAMMPS, although hopefully I can soon get that part working. It works on a single node environment, but as soon as you use MPI to add more cores we seg fault due to the LAMMPS stuff. Dont know why.


##Simulation checks
It is important to check on things as the simulation progesses. You should look at the V(k) for the most recent model file. The V(k) may be printed along with the model, but if not just run the model through femsim. You can change this in the code if you want (rmc.f90). 

Taking a look at the g(r) for a model is a great thing to do too.


##Autoslice
Autoslice, the slow but much more accurate intensity calculation, is implemented but I have no clue if it works correctly. Hopefully the default settings that are in there now are fine.

Femsim should always use autoslice; a simple if statement is all thats needed for this once its working.


##To list / notes
It may be necessary to reprogram the parallelization to work on any cluster. To do this I could start MPI on all the cores, but only do the main segment on one per processor. The other ones would just sit idle until they receive the "end" signal from the master process. Also, Odie has 8 cores per node, but thats in two quad-core processors. So OpenMP with 4 threads might be faster per core than OpenMP with 8 threads. The 4 cores per processor share fast cache RAM.

In the chi_squared subroutine do we (or should we) take into account data that is more difficult to obtain (such as that at lower k / q). I know we weight each type of data, but do we weight the data ranges themselves? Use Bayesian statistics to figure out the correct general weightings. The other approach is to use the refinement against only half the data. Then you simulate the other half and see whether it agrees - cross validation. Really really expensive. So do the first one - figure out how to and then do it.

Femsim should always use autoslice; a simple if statement is all thats needed for this once its working.
