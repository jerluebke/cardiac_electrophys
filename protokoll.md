# Content
1. Introduction
2. Three models
    2.1 Hudgkin & Huxley (1952)
    2.2 Aliev & Panfilov (1996)
    2.3 Fenton et al (2002)
3. Spatial dynamics in 1D
4. Spatial dynamics in 2D
5. Discussion
A. The code


# 1. Introduction
The human heart is a fascinating apparatus, which does its work in a
constant and reliable fashion - usually without disruption - for the whole
of a persons life. To put things into perspective, the heart of an average
human being performs
        70 (typical rest pulse per minute)
    x 1440 (minutes per day)
    x  365 (days per year)
    x   80 (estimated average human life time)
    ~ 3e9 beats per life time,
while a typical car engine performs
    300,000 (km driven during the cars life time)
    / 50    (km/h, average speed)   <- cars life time in hours
    x 60    (minutes per hour)
    x 2,200 (revolutions per minute)
    ~ 8e8 duty cycles per life time.

Investigating the hearts physical working principle poses an interesting
challenge with undoubtedly many relevant applications such as
gaining a deeper understanding of and developing more advanced treatments
for - many times dangerous - arrhythmia.

The heart muscles basic functionality is rhythmically contracting itself
triggered by electrical signals, which are being conducted by the heart
muscle cells (the cardiomyocytes) themselves (this is the remarkable thing
here, because usually this task is being taken care of by neural cells),
which means that this kind of tissue combines the ability to both perform
mechanical work and conduct electrical signals.

Going a little more into detail: pacemaker cells (specialized
cardiomyocytes) in the SA-node rhythmically generate action potential
which travel at about .05-1 m/s to the AV-node and from there after a
delay at about 2-4 m/s through the ventricular bundles and the Purkinje
fibres. Ultimately those action potentials cause the contraction of the
regions of heart tissue controlled by the respective bundles of conduction
cells.

## Cell structure and membrane potential
For a better understanding it is helpful to take a closer look at the
microscopic structure of the cardiomyocytes:
 * tubular cells containing chains of _myofibril_ (fibres composed of long
   proteins), which are responsible for contraction of the muscle tissue
 * _sarcoplastic reticulum_: membrane-enclosed regions, mainly storing 
   Ca^2+ ions
 * enclosed by a double lipid-layered membrane: the _sarcolemma_
In longitudinal direction _interlacing disks_ join the cells together and
via the _gap junctions_ allow propagation of action potentials. Because of
these features the heart muscle forms a _syncytium_, i.e. the single cells
behave like a single coordinated unit.

Now the interior and exterior (i.e. the intermediate space between
neighbouring cells) regions of a cells exhibit different concentrations of
various ion species (this imbalance is being maintained by special ion
pumps and gates in the cell membrane), which results in a voltage between
those regions: the membrane potential V=Φ\_i-Φ\_e, which in the rest case 
is equal to the rest potential V\_rest.

> TODO: table with typical values for V\_rest




>  vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : 
