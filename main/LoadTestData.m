function LoadTestData
%for now this just loads the data from plant.data, but later we can have
%this load a seperate file with some initial start time
global Plant TestData DateSim 
DateSim = Plant.Data.Timestamp(1);
TestData = Plant.Data;


