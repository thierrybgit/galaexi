# GPU Kernel Development for GALAEXI

This guide provides a step-by-step example of how to port a module, function or feature to use GPUs using the new multi-backend approach.
For illustrations, the ApplyJacobian module is used.

At the end of the guide is a reference for the device memory and state management APIs in GALAEXI, as well as references for GPU related 
data locations (module, include files, preprocessor directives, etc.)

```{admonition} How To Use This Guide
It is highly recommended that you **READ THROUGH THE ENTIRE GUIDE FIRST** before returning to the beginning and following the steps for your module/feature.
Reading through the whole process will allow you to see the whole picture before you start. This will help you make decisions in the early steps that will make the later steps easier.
```

## STEPS
### Step 1 -- Refactor existing code
As is described elsewhere in this Developer Guide, creating unit tests for new features is required in **FLEXI** and **GALÆXI**.
In this case, porting a feature counts as "creating" a new feature. If a unit test for the feature you are porting does not already exist,
you **MUST** create one.

```{admonition} Unit Tests For New Kernels
It is the policy of the the admins that any merge request submitted for a new ported feature in **GALÆXI** that does not
have an associated unit test will be rejected and remain so until a unit test is added. Zero exceptions. Period.
```

How to add a unit test is covered in detail in the chapter on [Testing](#09_testing).

Once you have a unit test for the feature you would like to port (whether it is new or pre-exisiting), a good starting point is to build and run the
code now, as is. This will confirm that you are starting from a working code and that the feature works before porting. This is important for sanity
checks and debugging later.

### Step 2 -- Refactor existing code
Now that you have a unit test for the feature that you are going to port to GPUs, you can actually start porting. This guide will use the *ApplyJacobian*
module as an example as it is simple and highlights almost all of the important features of the porting process.

At the beginning of our porting process the `applyjacobian.f90` and `applyjacobian.t90` files look like this:
```{figure} figures/12_kerneldev/applyjac_start.png
:align: center
:width: 85%
:name: applyjac_start

Starting point for the *ApplyJacobian* source code. `applyjacobian.f90` (left)  and `applyjacobian.t90` (right)
```

(_NOTE_: Clicking on the code pictures in this guide will enlarge them.)

The first thing to notice is that the *ApplyJacobian* module uses a pattern this guide will call "t90 templating". This pattern is common in **FLEXI** and uses two source codes files with the same name,
one with a `.t90` extension, the other with the `.f90` extension. Basically, the .f90 file includes a number of duplicate modules that each include the `.t90` file with a `#include` macro. t90 templating is an
unnecessary anti-pattern and should removed when encountered.

So the next step in your porting journey will be to remove this templating approach. If the module you are porting does not use t90 templating, you can skip this process.

```{admonition} Removing t90 Templating

To remove the t90 templating from the *ApplyJacobian* module, you will: 
1) Delete all of the modules in `applyjacobian.f90` except `MODULE MOD_ApplyJacobian`
2) Remove all instances of the `WITHnVar` preprocessor macros in both of our source code files. This functionally means that you will *always* pass the `TP_nVar` variable to the methods in the ApplyJacobian module and the multiple modules in `applyjacobian.f90` become redundant. For cleaner code, the function declaration lines in `applyjacobian.t90` can now be collapsed onto a single line.
3) Delete the line `#include "applyjacobian.t90` from `MODULE MOD_ApplyJacobian`
4) Paste the code from `applyjacobian.t90` into directly where the `#include "applyjacobian.t90` was. Make sure not to copy/paste the header comment or the `#include "flexi.h". They are already present in the other file.
5) Delete the file `applyjacobian.t90` as it is now completely empty and no longer needed

This process is the same for any other modules that use t90 templating.
```

Before proceeding any further you can examine the code further to see if more simplification is possible. **This highlights the fact that porting offers an opportunity to update old code.** ApplyJacobian is the perfect example of this. Looking through the three functions you now have in your `MOD_ApplyJacobian` module, you can see that these functions might be redundant. First off, a simple search of the code base for all calls to `ApplyJacobian` will reveal that the base overload, called just `ApplyJacobian`, is never invoked. So you can remove it. Easy. Next, you notice that `ApplyJacobian_select` is functionally identical to `ApplyJacobian_local` except
for a few additional lines of logic for FV sub-cells. This means that we can just remove `ApplyJacobian_local`. You now only have one method where you used to have three. This single method can do everything the other two functions could and even manages to still adhere to the "one function, one task" philosophy. Just to make your code that little bit simpler, you can even delete the `INTERFACE ApplyJacobian` at the top of the module and rename `ApplyJacobian_select` to `ApplyJacobian`. 

After all of that, your single source code file, `applyjacobian.f90`, should look like this:
```{figure} figures/12_kerneldev/applyjac_endstep1.png
:align: center
:width: 60%
:name: applyjac_endstep1

`applyjacobian.f90` after initial refactoring
```

```{admonition} Propogating Changes After Refactoring
*ApplyJacobian* is a core module and is called many times throughout the code. Since you have changed its API, there are now a number of calls to this subroutine that are no longer valid. You will have to go through the code and update those calls to reflect the new API.

Depending on which module and/or feature you are porting, this process may or may not be necessary.
```

Now that you have made signficant changes to the code, it would be a good idea to try to build and fix anything that breaks. You should also run your unit test to make sure the refactored method still produces the expected behavior.

At this point you might be thinking, "Where's the GPU? You said there would be GPU. You're a liar!" If that's the case, don't worry, you're getting there. These first two steps (unit tests and refactoring) are first for reason. You now have a way to test the important developments ahead and have confidence that they are correct when you finish. You also have made your life easier. By simplifying your module/feature, your code is easier to follow, easier to maintain and, in the case the *ApplyJacobian* module, now only one kernel is needed instead of three!

Unfortunately, before you do anything for GPUs, there are a few more steps that involve refactoring Fortran CPU code.

### Step 3 -- Create Fortran entry point method
**GALÆXI** supports *ALL* architectures, including CPUs. So since you already have a CPU backend with your existing Fortran function, the next step on your porting journey is to turn that method into a CPU backend.

```{admonition} Governing Philosophy of Porting in **GALÆXI**
Now is a good time to describe the core idea behind the decisions made in the design of **GALÆXI**'s multi-backend approach. This idea is simple, the barrier between the CPU and GPU is placed at the language barrier.

**From now on, if you are working on Fortran code, you are working on code that concerns the CPU. If you are working on C code, you are working on code that concerns the GPU.**

The key word in that manifesto is "concerns". Not all C code you see will *run* on the GPU. Sometimes C code will run on the CPU, but will always do something *concerning* the GPU, such as mutate or poll its current state (e.g. allocating memory on the GPU or checking to see how much GPU memory is available).
```

Returning to the example of the *ApplyJacobian* module, the first step in turn the existing Fortran method into a CPU backend is to introduce a new method into your module. This method will be called the "entry point" method from now on. This method's sole job is to prepare anything that is the same between the CPU and GPU backends and then launch the desired backend.

Again, this method will be added to the *ApplyJacobian* module and will look like this:
```{figure} figures/12_kerneldev/applyjac_entrypoint1.png
:align: center
:width: 60%
:name: applyjac_entrypoint1

"entry point" method example for *ApplyJacobian* backends
```
There are few features to recognize about this new "entry point" method:
1) It has the exact same name and arguments and uses the same modules as the `ApplyJacobian` method you already have.
2) It introduces an important preprocessor variable, `USE_ACCEL`. This variable is set by CMake during configuration and its use in the *ApplyJacobian* entry point will be a common pattern observed throughout your porting process. This variable will often be used in one of two ways:
   1) `#if (USE_ACCEL == ACCEL_OFF) ... #endif` -- code guarded by this preprocessor directive is only compiled when the **CPU** backends have been selected.
   2) `#if (USE_ACCEL != ACCEL_OFF) ... #endif` -- code guarded by this preprocessor directive is only compiled when the **GPU** backends have been selected.
3) Inside of the `#if (USE_ACCEL == ACCEL_OFF)` a function named `ApplyJacobian_Host` is called. From now on the `_Host` tag on a method will denote that that method is a CPU backend (The CPU is often called the "host" in GPU computing)

With this new entry point method added in its initial form, your code now has the machinery it needs to split between CPU or GPU backends, depending on how the code was compiled.

### Step 4 -- Refactor existing Fortran function as a CPU backend
Now looking at your entry point method, you will notice that the method `ApplyJacobian_Host` does not exist yet. Time to change that.

Going down to the "original" `ApplyJacobian` method (the one you created during the refactoring of Step 2), its time to go through the steps to transform this method into `ApplyJacobian_Host`. During this process its important to remember two rules for CPU backends in **GALÆXI**. These rules must be followed, no matter what.

```{admonition} Rules for CPU backends in **GALÆXI**
**RULE #1:** Using modules inside of CPU backends is forbidden. That is, having any occurrence of `USE MOD_???` is not allowed. This means that **_everything must be passed into the CPU backend as an argument_**. This will ultimately lead to some ugly function calls. That's fine. The goal of this rule is to make sure the function interfaces between ALL backends are the same. This makes debugging, maintaining and extending to other backends in the future easier.

**RULE #2:** At this point, the time for refactoring the Fortran code inside the CPU backend is over. That was finished in Step 2. DO NOT change the logic for the actual computation at this stage. You risk breaking something and then you won't be certain if what is wrong is your new interface or the changed logic. The process of software development is like performing an experiment, only change one independent variable at a time!
```

Following the above two rules, you can now refactor the "original" `ApplyJacobian` method into `ApplyJacobian_Host`, at the end of this process, it should look like this:

```{figure} figures/12_kerneldev/applyjac_hostbackend.png
:align: center
:width: 60%
:name: applyjac_hostbackend

Exisiting `ApplyJacobian` method refactored into a CPU backend
```

Again, notice that there are no `USE MOD_???` statements. Everything is passed as an argument now. You will also notice that the logic of the method is exactly the same as before.

And after you have refactored the backend method, you can update the call to it in the entry point.


### Step 5 -- Test CPU backend
You now have a CPU backend and the GPU stuff is just around the corner. Before you go onto to that though, now is time again to build, fix and test your new CPU backend. You could (and probably should) go beyond just running your unit test and see if you can run a test simulation with it. Before writing a single line of code concerning the GPU, everything concerning the CPU should be working perfectly. This also includes the new entry point method, which will be important for the GPU backend too.


### Step 6 -- Add interfaces to GPU backend on Fortran side
After the last step you now have a CPU-only code, just like you started with, but now it has been refactored such that it can accept a GPU backend.
So finally, its time.

What you came here for.

GPUs.

The first step in adding a GPU backend to this *ApplyJacobian* module is to create an interface for the GPU backend. This will be added to the top portion (before the `CONTAINS` statement) of the *ApplyJacobian* module and look like this:
```{figure} figures/12_kerneldev/applyjac_fortdevinterface.png
:align: center
:width: 60%
:name: applyjac_fortdevinterface

Exisiting `ApplyJacobian` Fortran interface to GPU backend
```
The important aspects of this new addition are:
1) Notice the pattern `#if (USE_ACCEL != ACCEL_OFF) ... #endif` that was mentioned earlier. This interface will only be compiled if you have turned on GPU acceleration during configuration. You must do this because the C code this interface will eventually link is not compiled when GPU acceleration is turned off. If you don't guard this interface with `#if (USE_ACCEL != ACCEL_OFF) ... #endif`, there will be errors in the linking stage when building **GALÆXI**.
2) The name of the function uses the `_Device` tag on the end, which is the GPU counterpart of the `_Host` tag seen earlier (the GPU is commonly referred to as the "device" in GPU computing).
3) The parameters passed to `ApplyJacobian_Device` are exactly this same as those passed to `ApplyJacobian_Host`. This is because the aforementioned desire to have the interfaces for all backends be consistent. However, you will notice there are two small differences with the parameters. The first is `streamID`. This is a new parameter only needed for the GPU side. It will be discussed more later in Step 8. The second is that `U` and `sJ` have become `d_U` and `d_sJ` and instead of `REAL` arrays, they are now `INTEGER` scalar values. The reasons for that will be explained shortly.
4) After the end of the `subroutine` declaration there is a new piece of code `BIND(C, NAME="ApplyJacobian_Device")`. This code points the compiler to the C function that this interface connects to and that we will write in Step 8. The value for `NAME` here should *always be exactly the same as the name of the interface it's connected to.* That isn't a requirement of the compiler, its just a requirement for code readability.
5) The interface includes a module `USE ISO_C_BINDING`. This module brings in some utilities needed for connected Fortran to C code that you need for this interface to function properly.
6) You will notice that, instead of plain `INTEGER`s like for the CPU backend, all of the values being passed into the GPU backend are of kind `INTEGER(C_INT)`. Adding in the `C_INT` isn't technically necessary, but it's good practice. You'll also notice an additional descriptor for those integers, `VALUE`. Fortran passes everything by reference by default (i.e. it passes pointers). Passing pointers between Fortran and C can (and usually does) cause a lot of headaches. Seeing as all that is happening here is the passing of a handful of single integer values, it is much easier (and perfectly acceptable) to pass by value and avoid the issue of passing pointers altogether.

Now that you have defined the interface to the GPU backend, you can add a call to it in the `#else /* USE_ACCEL */` branch of the entry point method `ApplyJacobian`.

But you cannot call this method yet. First off, the C method the interface links to doesn't exist yet. Second, the code doesn't know what the variables `d_U` and `d_sJ` are.

You will address the variables issue first.

```{admonition} Terminology
Going forward, the terms "CPU" and "GPU" will no longer be used. Instead, the terms "host" (for CPU) and "device" (for GPU) will be used. This is to be consistent with the terminology used in the **GALÆXI** code base (and literature on GPU accelerated computing)
```

### Step 7 -- Adding device variables and managing device memory
In the host backend, `d_U` and `d_sJ` were double precision real arrays, but in the device backend they are single integer values. Why this is brings up how device memory is handled in **GALÆXI**.

```{admonition} Device variables in **GALÆXI**
Device variables belong to the device. So following the core philosophy of **GALÆXI** from earlier, those variables must be defined and live on the C side of the code. This fact presents an issue, as there needs to be some way for the Fortran side (host) to tell the C side (device) which variables to use for each kernel. This problem is solved in **GALÆXI** with the use of a hash map data structure called `DeviceVars`. You will get to know this hash map well in the coming steps. When you allocate device memory, a pointer to that memory is stored into the `DeviceVars` hash map and the key to retrieve that pointer is passed back the Fortran side (host).

For consistency, ease and readability, all "device keys" (as they will be called from now on) should be named with the convention `d_<name_of_host_var>`. An example of this can be seeN in the *ApplyJacobian* module. There is a host-owned array called `U`. There also exists a device-owned copy of this array. The pointer to that device-owned copy is stored in the `DeviceVars` hash map. To key to fetch that pointer is saved in an integer variable called `d_U`.
```

The question from the above highlight is: Where does the device-owned copy come from and how is it associated to the variable `d_U`?

Enter the **GALÆXI** device memory management API! A full definition of this API is found [below](#MOD_DeviceMem).

To summarize, there exists a module in **GALÆXI** called `MOD_DeviceMem`. This module has a series of Fortran functions that can be used to allocate memory on the device and copy data back and forth between the host and device. These methods should be used in `Init` functions at the beginning of the simulation to allocate and initialize the device-owned copies of arrays.
For example, look how this is handled for the Jacobian array `sJ`:

```{figure} figures/12_kerneldev/sJ_devicemem.png
:align: center
:width: 60%
:name: sJ_devicemem

Allocation of device-owned copy of Jacobian array `sJ` and copy of initial values from the host to the device.
```

The above code exists in the `InitMesh` method found in `mesh.f90`. The call to `AllocateDeviceMemory` allocates a block of memory on the device with the same size as the host-owned array `sJ`. The pointer to the new GPU memory block is stored in the `DeviceVars` hash map and the key to retrieve it is stored in the `d_sJ` variable. You should always declare the device key variable along side its host counterpart (e.g. declare `d_sJ` along side `sJ` in `mesh_vars.f90`).

```{admonition} Note on using the device memory management API
All calls to the **GALÆXI** device memory managment API should be guarded behind `#if (USE_ACCEL != ACCEL_OFF) ... #endif` preprocessor flags. The device memory management API is not compiled for CPU builds of **GALÆXI**.
```

At this point in the development of **GALÆXI**, almost all core variables already have this allocate/copy process added. However, you should always check to make sure that the variable has had a device copy allocated and intialized. If it hasn't, add it now.

For sake of example, assume that there isn't a GPU-owned copy of `sJ` in the code already. You should now perform the following steps:
1) Add a declaration in the form `INTEGER(C_INT) :: d_sJ` to the module `MOD_Mesh_Vars` in `mesh_vars.f90`.
2) Add an allocation and copy step for the device-owned copies as illustrated in Fig. 6.

You can now add a guarded `USE` statement for `d_sJ` to the entry point method `ApplyJacobian`. For `d_U`, you will want to apply the Jacobian matrix to many different arrays throughout a computation, not just the conservative variables array (which is `U` represents). So for `d_U`, add an additional argument to the entry point method, which will allow calls to the entry point to pass keys for arbitrary arrays. Your entry point will now look like this:

```{figure} figures/12_kerneldev/applyjac_entrypoint_final.png
:align: center
:width: 60%
:name: applyjac_entrypoint_final

Finished form of the *ApplyJacobian* Fortran entry point method
```

You will notice there are few additional changes made to the code in Fig. 7. First, `streamID` is being passed as an argument to the entry point. This will allow different calls to `ApplyJacobian` to be run with different levels of priority on the GPU, a good feature to have. The other change is the addition of the line `USE MOD_Device`. That module contains a number of variables concerning host control of the device. In this case, you need access to the `STREAM_DEFAULT` variable to set a default priority level to `streamID`. A fair question to ask to that is, "What the hell is a stream? What are you talking about with 'priority levels'?". Unfortunately, that is a more advanced topic that won't be explained here. For now, just add a `streamID` variable to your entry point, set it to `STREAM_DEFAULT`, pass it to your device backend interface and forget about it. If you really *have* to know how streams work on GPUs, that's what the internet was invented for.

### Step 8 -- Add interfaces to GPU backend on C side
With the device backend interface added and the entry point method in a finished state, you are pretty much done with the Fortran side of this endeavor. Now it's time to start writing some C.

**NOTE** This guide won't even think about teaching you how write C code. That is the topic of literally entire stacks of books and cannot be achieved here. If you don't know how to write C, you should go do that now.

First thing you will do is create a new source code file for this C code. It should share the name of Fortran file all of the host code lives in and have the `.cu` extension. This new file should live along side its host counterpart in the same source directory. So for the *ApplyJacobian* module, you will create a new file called `applyjacobian.cu` in the same directory as `applyjacobian.f90`.

Before adding any code to the new source code file, go ahead and add the copyright header. This can just be copy and pasted from an existing `.cu` file elsewhere in the code.

The next thing you must do is create a C method for the `ApplyJacobian_Device` interface from Step 6 to attach to. Because you set `NAME="ApplyJacobian_Device` in `BIND(C)`, you must call this new C function `ApplyJacobian_Device`. The typing and names for the arguments to this new C function must match those of the Fortran interface.

```{admonition} Note on GPU programming terminology
Going forward in this guide, you will start to see the term "kernel" used quite a bit. A "kernel" is the code that actually runs on the device itself.
```

This new C function should look like this:
```{figure} figures/12_kerneldev/applyjac_kernel_launcher.png
:align: center
:width: 60%
:name: applyjac_kernel_launcher

Kernel launcher method for the *ApplyJacobian* module
```

A few things to point out about our new C method:
1) Even though this function is written in C, it actually runs on the host. Remember earlier it was mentioned that this would happen. Even though this function runs on the host, it *concerns* the device (in this case it changes the device's state by launching a kernel) and is therefore written in C.
2) The variable `nDOFs`. This variable is the number of degrees of freedom that the device kernel will perform calculations on. **GALÆXI**, in general, parallelizes kernels by DOF. This means that each thread on the device is given a single DOF to work on. You can see in the `INVOKE_KERNEL` macro that this `nDOFS` value is divided by 256 when passed. This is used in conjunction with the second passed value of 256 to determine the number of threads to use in the kernel. Because *ApplyJacobian* operates on the entire volume of the domain, `nDOFs` is set to the total number of DOFs in the domain.
3) The `INVOKE_KERNEL` macro. This macro takes formats its arguments into a device kernel call at compile time. This is done because different GPU programming models use different syntax for calling GPU code. These differences are hidden behind `INVOKE_KERNEL`. A full description of `INVOKE_KERNEL` is in the [API Reference section](#api-reference)
4) You can see in the arguments to `INVOKE_KERNEL` the device variable keys passed from the Fortran side, `d_U` and `d_sJ`, are used to index the `DeviceVars` hash map, just like was described in Step 7. Something additional here is the `(double *)` decorators added in front of `DeviceVars`. If you DON'T know C very well, just know that for double precision arrays like `U` and `sJ`, you have to include this decorator. If your Fortran array was an integer array (like is seen for `FV_elems` in the example), then you use `(int *)` instead. If you DO know C, then you can know that the device variable pointers are stored as `void *` in `DeviceVars` so only one hash map is needed. When those pointers are fetched, you have to cast the `void *` into the proper type.

From the above notes on the new C method, you will notice that all this method does is "invoke the kernel.", which is a fancy way of saying "it calls device code". From now on, this C method will be called the "kernel launcher".

Now that you have written the kernel launcher, you'd think it would be time to finally write the device code, but you'd be wrong! There is one more thing to do before that.

You see, you've been lied to. The code you've been writing hasn't actually been C code. It's been C++ code. Well, actually its been CUDA C++ code, but just pretend you didn't read that. The reason it has been called C code so far is because, functionally, it is. Many of the fancy features of C++ (like classes and the C++ STL) are not supported in kernel code, so it's best to just think about it as C code. It makes things easier in the long run.

Regardless of how you think about your code, it still is *technically* C++ and because of that you must include what will be called the "extern C prototype". This is just a function prototype for the kernel launcher wrapped in `extern "C"{}`. The reasons for this are technical and not necessary to know at this point. Just do it.

The extern C prototype should look like this:
```{figure} figures/12_kerneldev/applyjac_externc.png
:align: center
:width: 60%
:name: applyjac_externc

extern C prototype for ApplyJacobian module
```

```{admonition} Order of functions in the C code
The extern C prototype for the kernel launcher must be placed at the very top of the `.cu` file (directly below the copyright header). The kernel launcher must be placed at the very bottom. These are requirements for the compiler.
```

Just in case you are lost or confused as to what all of these launchers, interfaces, prototypes, etc. are doing, here is a quick overview of the code flow for the device backend:

1) Some function in the main Fortran code calls the `ApplyJacobian` entry point on the host
2) The entry point calls `ApplyJacobian_Device`
3) The call to `ApplyJacobian_Device` routes through the Fortran interface of the same name, through the extern C prototype to the kernel launcher.
4) The kernel launcher launches the kernel `ApplyJacobian_Kernel` on the device.

The only missing link in this chain now is the kernel itself.


### Step 9 -- Write kernel
Explaining this is easier if you first look at the example kernel below:
```{figure} figures/12_kerneldev/applyjac_kernel.png
:align: center
:width: 60%
:name: applyjac_kernel

Example device kernel for the *ApplyJacobian* module
```

**NOTE** The logic for handling FV sub-cells was removed from Fig. 10 for simplicity

Here is everything you should take note of from the kernel code:
1) The computation logic in this "initial version" of the kernel flows exactly like the code in the Fortran-based host backend `ApplyJacobian_Host`. This is pointed out to emphasize that, when writing the C kernel for the first time, try to follow the original Fortran as closely as possible. What results might not be "optimum" C code, but it should work with only a few changes and tweaks. This code can be updated and optimized later on to make it more "C like".
2) The function has the `__global__` decorator. This decorator tells the compiler that the kernel function is a function that runs on the device but is called from the host. Throughout existing kernels in **GALÆXI**, there are also some examples of functions that use the `__device__` decorator. Those functions run on the device and are also called from other device functions (that is, they are completely unreachable from the host).
3) The line `int threadID = blockIdx.x * blockDim.x + threadIdx.x;`. This line uses the CUDA/HIP API to determine the global index of the thread currently running the function. How exactly this line works isn't important really, just copy and paste into kernels where needed and know that it returns what it returns. Because you launched the kernel with a single thread assigned to a single DOF, this `threadID` is also the global ID of the DOF the thread will work on.
4) The line `if (threadID < nDOFs)`. This if statement filters out any device threads that have been assigned a DOF that doesn't exist. This can happen because the device launches threads in numbers that are typically powers of 2. If you have a mesh that has a number of DOFs that isn't nicely divided by powers of 2, it won't be possible to launch the exact number of threads needed for the number of DOFs you have. To overcome this, more threads than are needed are initially launched and then those assigned DOFs that are "out of bounds" simply aren't allowed to perform any computation.
5) The first line inside of the if statement is `U_offset = threadID*TP_nVar;`. This brings up an incredibly important fact about how memory is handled inside kernels in **GALÆXI**. **_Multi-dimensional arrays in the C code are treated as flattened 1D arrays_**. Fortran arrays, regardless of what dimensionality they are assigned, are stored in memory as a single long line of bytes. The order of these bytes is determined by the Fortran standard. It is actually much, much easier for device code to use knowledge of how multi-dimensional Fortran arrays are stored in memory and index them as 1D memory. Interfacing Fotran arrays to multi-dimensional arrays in C can be a tremendous pain (and C does not handle multi-dimensional arrays well anyway). Indexing flattened memory can be a difficult thing to wrap one's mind around at first, but you will get a feel for it. Luckily, there are number of functions already available in **GALÆXI** to calculate the offsets in flattened memory for you. They are described in the [API reference section](#offsetting-methods-for-flattened-memory).
6) Another thing that can be very difficult to wrap one's mind around is this continuing idea that one DOF is assigned to one device thread. You will notice that there are no `for` loops in the kernel code. That is because we don't need them. The kernel code is written as if it only handles a single DOF. Think of a kernel as the code from *inside* a `for` loop and the GPU is looping through all of the DOFs for us, running the kernel function `nDOFs` times, once for each DOF. Understanding this idea is the crux of programming GPUs. Once you get that concept, this will all get much easier.

To finish up the kernel writing process, there is one last thing to do. There is a file in `$GALAEXI_ROOT/src/device` called `main.cu`. In that file is a list of `#include` statements with the names of `.cu` files. If you created a new `.cu` file, add a new `#include` line for it in that file. Why this is done has to do with how the CUDA and HIP compilers work and how they optimize and link code. At this point, you only need to know that this process makes the code faster. When placing the new `#include` in `main.cu` please note the order of the `#include` lines and place your new file under any other modules that it might depend on (this file is NOT dependency checked).

### Step 10 -- Test device backend
So you've just written a kernel. Awesome. Before getting too ahead of yourself though, you need to do a 5 more things:
1) Update that unit test from Step 1 to test your device backend.
2) Run that unit test on your kernel and fix everything that breaks (there will probably be a lot. There always is).
3) Rerun a simple test problem with GPU acceleration active to make sure everything is performing as expected.
4) If you developed the new kernel on an NVIDIA-based system, you should then go test your kernel on an AMD-based system.
5) Rebuild and retest the CPU backend again, just to make sure nothing regressed.

```{admonition} Hooray!
Alright, congratulations. You're done. You have ported a feature in **GALÆXI**!

Push it to the GitLab and save time later to completely redo it all after the reviewer of your merge request rejects it!
```
---
## API reference
### device.h -- C side device variables and macros definitions
The following variables and macros are found in `$GALAEXI_ROOT/src/device/device.h`. This header file is intended for use **ONLY IN C FILES**.

`$GALAEXI_ROOT/src/device/device.h` is included in `main.cu`, making these variables macros available to any `.cu` file included in `main.cu`.

```{eval-rst}
.. cpp:struct:: DeviceMemoryState

   Simple struct that stores return values of a call to (cuda/hip)GetMemInfo

   .. cpp:member:: size_t FreeMemory

      Amount of free memory currently available on the device

   .. cpp:member:: size_t TotalMemory

      Amount of total available memory on the device


.. cpp:var:: std::unordered_map<int, void*> DeviceVars

|  This hash map stores pointers to all variables in device memory. It is keyed by and integer value stored in a variable with the name of the device variable with the pattern
|  host variable       --> var
|  device variable key --> d_var
|  Instantiated in device_mem.cpp. Only has scope within libflexif90 compile unit (can't be used in the entry point method from flexi.f90).


.. c:macro:: DEVICE_VAR_CHECK( call )

   Wrapper macro that automatically inserts error handling code around CUDA/HIP API calls.

   :param call: CUDA/HIP API function call.


.. c:macro:: INVOKE_KERNEL( func, blocks, threads, shared, stream, ... )

   Wrapper macro that formats kernel launches into the appropriate format for a given vendor. For example a CUDA kernel launch will be formatted as func<<<blocks,threads,shared,stream>>>(__VA_ARGS__). A semi-colon (;) must be added to the line that uses INVOKE_KERNEL, the macro doesn't add it.

   :param func: Name of the kernel function
   :param blocks: Number of thread blocks to launch the kernel with
   :param shared: Amount of shared memory to reserve for the kernel. Should almost always be passed as "0"
   :param stream: Stream ID to run the kernel on
   :param ...: After stream should follow the comma separated list of arguments for the kernel function itself
```
---
### Offsetting methods for flattened memory
The following methods are found in `$GALAEXI_ROOT/src/device/offsets.cu`. They are intended for finding offsets for flattened memory inside of kernels.

`$GALAEXI_ROOT/src/device/offsets.cu` is included in `main.cu`, making these methods available to any `.cu` file included in `main.cu`.

```{eval-rst}

.. cpp:function:: int FindElemOffset(int nVar, int Nloc, int elemIdx)

   Find the offset for a specific element in 1D flattened device mem

   :param nVar: Number of flow variables
   :param Nloc: Local polynomial order
   :param elemIdx: Index of the desired element (assumed to be indexed from 1)
   :returns: Computed offset in bytes to the starting index of the desired element

.. cpp:function:: int FindElemMetricsOffset(int Nloc, int elemIdx, int FV_idx)

   Find the offset of an element in the Metrics arrays

   :param Nloc: Local polynomial order
   :param elemIdx: Index of the current element (assumed to be indexed from 1)
   :param FV_idx: Idx for FV_SIZE  (assumed to be indexed from 0)

.. cpp:function:: int FindElemJacobianOffset(int Nloc, int elemIdx, int FV_idx)

   Find the offset of an element in the Jacobian array

   :param Nloc: Local polynomial order
   :param elemIdx: Index of the current element (assumed to be indexed from 1)
   :param FV_idx: Idx for FV_SIZE (assumed to be indexed from 0)

.. cpp:function:: int IndexFlatFortranArr(int xSize, int xIdx, int yIdx)

   Overload for 2D memory that is allocated from 1 in both dimensions. ASSUMES THAT xIdx AND yIdx ARE INDEXED STARTING AT 1.

   :param xSize: Size of the 1st dimension of the array
   :param xIdx: Fortran index of the 1st dimension of the array you would like the offset for (assumed to be indexed from 1)
   :param yIdx: Fortran index of the 2nd dimension of the array you would like the offset for (assumed to be indexed from 1)
   :returns: Calculated offset of the desired element in flattened memory

.. cpp:function:: int IndexFlatFortranArr(int xSize, int ySize, int xIdx, int yIdx, int zIdx)

   Overload for 3D memory that is allocated from 1 in all dimensions. ASSUMES THAT xIdx, yIdx AND zIdx ARE INDEXED STARTING AT 1.

   :param xSize: Size of the 1st dimension of the array
   :param ySize: Size of the 2nd dimension of the array
   :param xIdx: Fortran index of the 1st dimension of the array you would like the offset for (assumed to be indexed from 1)
   :param yIdx: Fortran index of the 2nd dimension of the array you would like the offset for (assumed to be indexed from 1)
   :param zIdx: Fortran index of the 3rd dimension of the array you would like the offset for (assumed to be indexed from 1)
   :returns: Calculated offset of the desired element in flattened memory

.. cpp:function:: int FindPointOffset(int nVar, int Nloc, int elemIdx, int i, int j, int k)

   Find the offset for a specific DOF in 1D flattened device mem

   :param nVar: Number of flow variables
   :param Nloc: Local polynomial order
   :param elemIdx: Index of the desired element (Indexed from 1)
   :param i: Index of the first dimension DOF of the element. (Indexed from 0)
   :param j: Index of the second dimension DOF of the element. (Indexed from 0)
   :param k: Index of the third dimension DOF of the element. (Indexed from 0)
   :returns: Computed offset in bytes to the starting index of desired single DOF

.. cpp:function:: int FindPointMasterSlaveOffset(int nVar, int Nloc, int sideIdx, int i, int j)

   Find the offset of a point in a flattened master/slave sides array

   :param nVar: Number of variables at each DOF
   :param Nloc: Local polynomial order
   :param sideIdx: Index of the desired side (assumed to be indexed from 0)
   :param i: Index of the first dimension DOF of the element. (assumed to be indexed from 0)
   :param j: Index of the second dimension DOF of the element. (assumed to be indexed from 0)
   :returns: Computed offset in bytes to the starting index of the desired DOF

.. cpp:function:: int FindPointMetricsOffset(int Nloc, int elemIdx, int nElems, int FV_idx, int i, int j, int k)

   Find the offset of a point in a flattened metrics array

   :param Nloc: Local polynomial order
   :param elemIdx: Index of the desired element (assumed to be indexed from 1)
   :param FV_idx: Idx for FV_SIZE  (assumed to be indexed from 0)
   :param i: Index of the first dimension DOF of the element. (assumed to be indexed from 0)
   :param j: Index of the second dimension DOF of the element. (assumed to be indexed from 0)
   :param k: Index of the third dimension DOF of the element. (assumed to be indexed from 0)
   :returns: Computed offset in bytes to the starting index of the desired DOF

.. cpp:function:: int FindPointJacobianOffset(int Nloc, int elemIdx, int nElems, int FV_idx, int i, int j, int k)

   Find the offset of a point in a flattened Jacobian array

   :param Nloc: Local polynomial order
   :param elemIdx: Index of the desired element (assumed to be indexed from 1)
   :param FV_idx: Idx for FV_SIZE (assumed to be indexed from 0)
   :param i: Index of the first dimension DOF of the element. (assumed to be indexed from 0)
   :param j: Index of the second dimension DOF of the element. (assumed to be indexed from 0)
   :param k: Index of the third dimension DOF of the element. (assumed to be indexed from 0)
   :returns: Computed offset in bytes to the starting index of the desired DOF

.. cpp:function:: int FindPointMappingOffset(int sideID, int Nloc, int i, int j, int k, int flip, bool isMapping2)

   Find the offset of a point in one flattened mapping arrays (S2V, V2S, S2V2, etc.)

   :param sideID: Index of the current side of the desired element. (assumed to be indexed from 0)
   :param Nloc: Local polynomial order 
   :param flip: Index of flip side (assumed to be indexed from 0)
   :param i: Index of the first dimension DOF of the element. (assumed to be indexed from 0)
   :param j: Index of the second dimension DOF of the element. (assumed to be indexed from 0)
   :param k: Index of the third dimension DOF of the element. (assumed to be indexed from 0)
   :param isMapping2: Is the mapping array being indexed S2V2 (true) or S2V (false)?

.. cpp:function:: int FindSideOffset(int nVar, int Nloc, int elemIdx, int sideIdx)

   Find the offset for a specific side of an element in 1D flattened device mem. This method assumes that the sideID is indexed from 0.

   :param nVar: Number of flow variables
   :param Nloc: Local polynomial order
   :param elemIdx: Index of the desired element  (assumed to be indexed from 1)
   :param sideIdx: Index of the side. (assumed to be indexed from 0)
   :returns: Computed offset in bytes to the starting index of the side within the desired element

.. cpp:function:: int FindSideMasterSlaveOffset(int nVar, int Nloc, int sideIdx)

   Find the offset to the start of a side in a master/slave sides array

   :param nVar: Number of variables at each DOF
   :param Nloc: Local polynomial order
   :param sideIdx: Index of the desired side (assumed to be indexed from 0)
   :returns: Computed offset in bytes to the starting index of the desired side

.. cpp:function:: int FindSideSurfaceDataOffset(int nVars, int Nloc, int nSides, int FV_idx, int sideIdx)

   Find the offset to the start of a side in the surface vector arrays

   :param nVars: Number of vars in the 1st dim in the array. Used to allow this method to offset arrays like NormVec (nDim=3) and SurfElem (nDim = 1)
   :param Nloc: Local polynomial order
   :param sideIdx: Index of the side to offset to (assumed to be indexed from 0)
   :returns: Computed offset in bytes to the starting index of the desired side

.. cpp:function:: ThreadIndicesSide FindThreadIndicesSide(int Nloc, int startSide)

   | Wraps the often used logic of finding the side and DOF indices for a device thread.
   | WARNING! Indexes threadIDs from 1, not 0. This is due to the logic for finding p and q.

   :returns: An instance of the ThreadIndices class initialized with the side ID and DOF indices of the calling device thread.
   :param Nloc: Local polynomial order
   :param startSide: For kernels that don't work on a full sides array, this allows offseting SideID.

.. cpp:function:: ThreadIndicesVolume FindThreadIndicesVolume(int Nloc)

   | Wraps the often used logic of finding the element and DOF indices for a device thread.
   | **WARNING!** Indexes threadIDs from 1, not 0. This is due to the logic for finding i,j,k.

   :returns: An instance of the ThreadIndices class initialized with the element ID and DOF indices of the calling device thread.
   :param Nloc: Local polynomial order

.. cpp:Struct:: ThreadIndicesSide

   Holds index information about the specific DOF assigned to a device thread in a kernel that operates on a single side. It hides the logic for calculating these IDs and indices with the method FindThreadIndicesSide.

   .. cpp:member:: int threadID

      ID of the current device thread. Indexed from 1.

   .. cpp:member:: int SideID

      ID of the side that has the DOF assigned to this thread. Indexed from 1.

   .. cpp:member:: int i, j, k

      Indices within side SideID for the DOF assigned to thread. Indexed from 0.

.. cpp:struct:: ThreadIndicesVolume

   Holds index information about the specific DOF assigned to a device thread in a kernel that operates on the whole volume. It hides the logic for calculating these IDs and indices with the method FindThreadIndicesVolume.

   .. cpp:member:: int threadID

      ID of the current device thread. Indexed from 1.

   .. cpp:member:: int ElemID

      ID of the element that has the DOF assigned to this thread. Indexed from 1.

   .. cpp:member:: int i, j, k

      Indices within element ElemID for the DOF assigned to thread. Indexed from 0.

```

---
### MOD_Device -- Resources for device handling on Fortran side

```{eval-rst}

**Module** :f:mod:`MOD_Device`

.. f:module:: MOD_Device
   :synopsis: This Fortran module holds copies or interfaces to C++ device variables in device.h needed on the Fortran side.

.. f:variable:: CurrentDeviceMemState
   :type: CheckDeviceMemoryState

   Instance of the DeviceMemoryState type that can be used globally to check state and then keep it stored, so it can be checked without having to repoll if we that state hasn't changed.

.. f:subroutine:: DefineParametersDevice()

   Define device input parameters

.. f:type:: DeviceMemoryState

   Mirrors the C struct of the same name defined in device.h. Only used as the return value for CheckDeviceMemoryState

   :f INTEGER(C_SIZE_T) FreeMemory: Amount of free memory currently available on the device

   :f INTEGER(C_SIZE_T) TotalMemory: Amount of total available memory on the device

.. f:subroutine:: InitDevice()

   Intialize the device

.. f:subroutine InitStreams(useStreams)

   Intializes device streams.

.. f:subroutine FinalizeDevice()

   Destroys any intialized streams at the end of the run.

.. f:variable:: STREAM_DEFAULT
   :type: INTEGER
   :attrs: PARAMETER = 0

.. f:variable:: STREAM_LOW
   :type: INTEGER
   :attrs: PARAMETER = 1

.. f:variable:: STREAM_MID
   :type: INTEGER
   :attrs: PARAMETER = 2

.. f:variable:: STREAM_HIGH
   :type: INTEGER
   :attrs: PARAMETER = 3

```

---
### MOD_DeviceManage -- Device management API

```{eval-rst}

**Module** :f:mod:`MOD_DeviceManage`

.. f:module:: MOD_DeviceManage
   :synopsis: Fortran interfaces to C++ methods used to issues commands to the device or change its operational state.

.. f:subroutine:: SetDevice(deviceID)

   Wrapper for the set device methods in CUDA and HIP. For now just used to create device contexts in unit tests.

   :param INTEGER(C_INT) deviceID [VALUE]: ID of the device being chosen

.. f:subroutine:: SynchronizeDevice()

   Wrapper for calls to CUDA/HIP API device synchronize methods. Before synchronizing the device, peeks for errors in recent kernel executions and throw exceptions if kernel errors are detected.

```
---
### MOD_DeviceMem -- Memory management API

```{eval-rst}

**Module** :f:mod:`MOD_DeviceMem`

.. f:module:: MOD_DeviceMem
   :synopsis: Fortran interfaces to C++ device memory methods

.. f:subroutine:: AllocateDeviceMemory(dVarKey, typeSize_bytes, arraySize)

   Allocate a block of memory on the device and store the pointer in the device vars map using the requested key.

   :param INTEGER(C_INT) dVarKey [VALUE]: Desired device variable key
   :param INTEGER(C_SIZE_T) typeSize_bytes [VALUE]: Size in bytes of the data type to be allocated
   :param INTEGER(C_INT) arraySize [VALUE]: Number of elements of type to allocate (1 for scalars)

.. f:subroutine:: CheckDeviceMemoryState(FreeMemory, TotalMemory)

   Check the current memory quota of the device

   :param INTEGER(C_SIZE_T) FreeMemory: Amount of free memory currently available on the device
   :param INTEGER(C_SIZE_T) TotalMemory: Amount of total available memory on the device

.. f:subroutine:: CopyFromDevice(hostMemory, dVarKey, arraySize)

   Copy data from the device to the host

   :param INTEGER(*)/DOUBLE(*) hostMemory: Host copy of the data
   :param INTEGER(C_INT) dVarKey [VALUE]: Device variable key
   :param INTEGER(C_INT) arraySize [VALUE]: Number of elements allocated for the block being copied (1 for scalars)

.. f:subroutine:: CopyToDevice(dVarKey, hostMemory, arraySize)

   Copy data to the device from the host

   :param INTEGER(C_INT) dVarKey [VALUE]: Device variable key
   :param INTEGER(*)/DOUBLE(*) hostMemory: Host copy of the data
   :param INTEGER(C_INT) arraySize [VALUE]: Number of elements allocated for the block being copied (1 for scalars)

.. f:subroutine:: FreeDeviceMemory(dVarKey)

   Deallocate memory on the device

   :param INTEGER(C_INT) dVarKey [VALUE]: Device variable key

.. f:variable:: SIZE_C_DOUBLE
   :type: INTEGER(C_SIZE_T)
   :attrs: PARAMETER = 8_C_SIZE_T

.. f:variable:: SIZE_C_INT
   :type: INTEGER(C_SIZE_T)
   :attrs: PARAMETER = 4_C_SIZE_T

```