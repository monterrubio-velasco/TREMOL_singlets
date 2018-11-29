# Author: Stefan Wiemer 03/95 

function  bmemag(VecMag)

# bmemag.jl: Function calculates the mean magnitute, the b value based on the mean and the standart deviation 
 
# Parameters
# - `VecMag ::Array`: Vector that contains the equivalent magnitudes

# return
# - `meanm1::Float`: mean value of VecMag vector
# - `b1::Float`: b-value 
# - `sig1::Float`: standard deviation
# - `av2::Float`: a-value

    newcat = copy(VecMag)
    maxmag = maximum(newcat)
    mima = minimum(newcat)
    if mima > 0
        mima = 0
    end
  
  #calculate the mean magnitude, b(mean) and std    
    n = length(newcat)
    meanm1 = mean(newcat)
    b1 = (1/(meanm1-minimum(newcat-0.05)))*log10(exp(1))
    sig1 = (sum((newcat-meanm1).^2))/(n*(n-1))
    sig1 = sqrt(sig1)
    sig1 = 2.30*sig1*b1^2            # standard deviation
    av2 = log10(length(newcat))+b1*minimum(newcat)
    return meanm1, b1, sig1, av2
end
