# MotionPlanningControl_git
Codes for paper "Bottom-level motion control for robotic fish to swim in groups: modeling and experiments"

Run "MotionPlaningControl.m" for Fig9, "MotionPlaningControl2.m" for Fig10, and "MotionPlaningControl3.m" for Fig11.

If you just want to use our model, you just need to use this part:
```MATLAB
    sol=ode45(@Equation2 ,[t ,t + dt] ,...
        [0.12, 0 ,0, 1.0, X,108680.61e-7]);
    ndt = size(sol.y,2);
    for ii = 1:size(sol.y,2)
        p.w(end+1)=p.w(end)+dt/ndt*deval(sol,sol.x(ii),3);
        p.x(end+1)=p.x(end)+dt/ndt*cos(p.w(end))*(deval(sol,sol.x(ii),1) - 0.12)...
            -dt/ndt*sin(p.w(end))*deval(sol,sol.x(ii),2);
        p.y(end+1)=p.y(end)+dt/ndt*sin(p.w(end))*(deval(sol,sol.x(ii),1) - 0.12)...
            +dt/ndt*cos(p.w(end))*deval(sol,sol.x(ii),2);
    end
```


## References:

If you use this code or data we kindly as that you please [cite Li et al, 2019](https://doi.org/10.1088/1748-3190/ab1052) 


Please check out the following references for more details:

    @article{li2019bottom,
            title={Bottom-level motion control for robotic fish to swim in groups: modeling and experiments},
            author={Li, Liang and Liu, Anquan and Wang, Wei and Ravi, Sridhar and Fu, Rubin and Yu, Junzhi and Xie, Guangming},
            journal={Bioinspiration \& biomimetics},
            volume={14},
            number={4},
            pages={046001},
            year={2019},
            publisher={IOP Publishing}}
