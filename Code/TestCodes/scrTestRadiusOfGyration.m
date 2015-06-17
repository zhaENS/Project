% scrTestRadiusOfGyration

% see that the chain occupies a region proportional to the radius of
% gyration sqrt(numBeads/6)*b

  % Initialize simulator framework parameters
    simulatorParams = SimulationFrameworkParams('showSimulation',true,...
                                                'numSteps',1000,...
                                                'dimension',3,...
                                                'dt',0.3,...
                                                'objectInteraction',false);      
                                            
                                            
    % % create a chain
    chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,...
                                     'springForce',true,...
                                     'bendingElasticityForce',false,...
                                     'bendingConst',1,...
                                     'springConst', 1,...
                                     'openningAngle',pi,...
                                     'minParticleEqDistance',0);
    
    cp          = ChainParams('numBeads',1000,...
                              'dimension',simulatorParams.simulator.dimension,...
                              'initializeInDomain',1,...
                              'forceParams',chainForces,...                              
                              'b',sqrt(3));
                                                                                        
    % create an spherical domain
    sphereForce = ForceManagerParams('lennardJonesForce',false,...
                                         'LJPotentialWidth',0.1,...
                                         'LJPotentialDepth',0.1,...
                                         'diffusionForce',true,...
                                         'diffusionConst',1,...
                                         'mechanicalForce',false,...
                                         'mechanicalForceDirection','out',...
                                         'mechanicalForceCenter',[0 0 0],...
                                         'mechanicalForceMagnitude',0,...
                                         'dt',simulatorParams.simulator.dt);
    
    dp(1)         = DomainHandlerParams('domainShape','sphere',...
                                        'reflectionType','off',...
                                        'dimension',simulatorParams.simulator.dimension,...
                                        'domainCenter',[0 0 0],...
                                        'forceParams',sphereForce,...                                        
                                        'domainWidth',sqrt(cp.numBeads/6)*cp.b,...
                                        'dimension',simulatorParams.simulator.dimension);
                                    
   % register the object parameters in the simulator framework
    simulatorParams.SetDomainParams(dp);
    simulatorParams.SetChainParams(cp);
    
    % Initialize simulator framework
    r = RouseSimulatorFramework(simulatorParams);
    % Run until relaxation time
    r.Run                                                                                                            