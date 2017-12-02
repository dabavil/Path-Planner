### This document describes the model used to plan and generate trajectories in this path planning project

## Result

The model was able to reliably drive the car around the track, exceeding 10 miles without accident:

![Screenshot of simulator](/screen_cap.png)


## Architecture

To deliver this result, the model goes throught the following steps, all of which are described in more detail below:

1. Get clarity on the situation with other cars (from line 248)
2. Evaluate the situational data (from line 297)
3. Calculate the cost of each option and decide (from line 364)
4. Generate trajectories


## More detailed description of the steps

# Get clarity on the situation with other cars 

In this section, we cycle through the information on all other cars in the vicinity, and exctracts how far ahead they are, how quickly they go, and in which lane they are. For each car this information flows further to step 2 (evaluation)


# Evaluate the situational data

In this section, the model uses the information about each car, to establish what is the situation in each lane. In particular: a) is it safe to move into this lane, or is there a risk of collision? and b) looking ahead, what is the max speed ego car will be able to go in this lane? If no cars are visible, the lane speed equals to car's max allowed speed (which is just under the 50 MPH limit)


# Calculate the cost of each option and decide

Here, the model evaluates a simple cost function, that takes into account both the information we know about the lane, and the current situation of the car. The cost increases if moving into the lane would present a collision danger (with a very large penalty). The cost also increases, if the speed possible in the lane is lower than the maximum speed. Finally, there is cost associated with each lane switch, to prevent the car from oscilating between the lanes too frequently. Costs were carefully tweaked to provide the result, and with simple tweaking of the speed and lane change penalties the behavior (or agressiveness) of the car can be altered dramatically. Once the evaluation is complete, the car moves one lane over in the direction of the optimal lane (unless already there)


# Generate trajectories

The core of this section is using Spline and related techniques revealed in the walkthrough video, to create a smooth, jerk-minimizing trajectory. I've tweaked several of the hyperparameters to achieve better performance. In particular - the distance of the seed points for spline (lines 481-583) are dynamically adjusted based on the current speed of the car. This allows the ego car to change lanes nimbly in slow traffic, which helps prevent collisions, yet results in smooth trajectories generated at higher speeds, thus preventing max jerk warnings. 

## Reflection

The amount of situations that the car can get into is quite surprising, even though the setup of the simulator is relatively simple, and the information available is near perfect (E.g., other cars ground truth). As a result, I found that while it's easy to make the car go the first mile, to acually complete a lap the model has to be relatively robust. My particular implementation is quite procedural, and has currently lots of potential to refactor and streamline the code (#TODO)




