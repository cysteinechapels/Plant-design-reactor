function OptimalFlame

LB = [ 0 0 0 0]
UB = [ 1 1 1 1]

Initial = [0.5 0.2 0.1 0.2]


function error = goal(x)
        [selectivity,O2perc,i]=Flammability(x);
        error = -selectivity;
    end


S = fmincon(@(x) goal(x),Initial,[],[],[],[],LB,UB);

S

[selectivity, O2perc, i] = Flammability(S)
end