import numpy
import matplotlib
import matplotlib.pyplot as plt; plt.rcdefaults()
matplotlib.rcParams['text.usetex'] = True

f2_scale = numpy.array([-17.9200220901249930755098639008252582465,
                        1.845764788512892377453602431656624264454,
                        1.839611388308253894657545680555773464748,
                        1.840275804402260393066356407589470103773,
                        1.840265160571370863263682954232296794854,
                        1.840265762727976760495004292302021649702,
                        1.84026576316694892619408592555747447912,
                        1.840265763131554575497970641550901728001
                        ])

Delta_constraint = [[0, 0.1, 1],
                    [0, 0.1, 1],
                    [0, 0.1, 1, 1.044339780356917847849369292480503379774],
                    [0, 0.1, 1, 1.044339780356917847849369292480503379774],
                    [0, 0.1, 1, 1.042500153473850784545700476498155026252,
                     1.044339780356917847849369292480503379774],
                    [0, 0.1, 1, 1.042500153473850784545700476498155026252,
                     1.044339780356917847849369292480503379774],
                    [0, 0.1, 1, 1.042500153473850784545700476498155026252,
                     1.044339780356917847849369292480503379774],
                    [0, 0.1, 1, 1.042496785729767303511021191519693984244,
                     1.042500153473850784545700476498155026252,
                     1.044339780356917847849369292480503379774]
                    ]

fig, ax = plt.subplots()
ax.set_xlim(0,1.5)
ax.set_ylim(-1,10)

Delta_0 = numpy.arange(0,2,0.01)
f1_0 = 1 + numpy.power(Delta_0,4)
f2_0 = numpy.power(Delta_0,2) + numpy.power(Delta_0,4)/12

ax.plot(Delta_0, f1_0 - f2_scale[0] * f2_0, linestyle="dotted", label="Iteration 1", color="xkcd:purple")
ax.plot(Delta_0, f1_0 - f2_scale[1] * f2_0, linestyle="--", label="Iteration 2", color="xkcd:green")
ax.plot(Delta_0, numpy.zeros(len(Delta_0)), color="xkcd:black")

ax.plot(numpy.array(Delta_constraint[0]),
        numpy.zeros(len(Delta_constraint[0])),"o",
        label="Initial Points", color="xkcd:blue")
ax.plot(Delta_constraint[2][3], 0, "s", label="First Added Point", color="xkcd:red")

ax.set_ylabel('$F_0$', fontsize=15)
ax.set_xlabel('$\Delta$', fontsize=15)
# ax.set_title('Time for Second Iteration',fontsize=18)

ax.legend()
plt.savefig("outer_plots_1.png",dpi=300)

plt.clf()
fig, ax = plt.subplots()
ax.set_xlim(1.025,1.06)
ax.set_ylim(-0.0001,0.0015)

Delta_1 = numpy.arange(1.02,1.06,0.00002)
f1_1 = 1 + numpy.power(Delta_1,4)
f2_1 = numpy.power(Delta_1,2) + numpy.power(Delta_1,4)/12

ax.plot(Delta_1, f1_1 - f2_scale[2] * f2_1, ":", label="Iteration 3", color="xkcd:brown")
ax.plot(Delta_1, f1_1 - f2_scale[3] * f2_1, linestyle=(0, (5, 1)), label="Iteration 4", color="xkcd:pink")
ax.plot(Delta_1, numpy.zeros(len(Delta_1)), color="xkcd:black")

ax.plot(Delta_constraint[2][3], 0, "s", label="First Added Point", color="xkcd:red")
ax.plot(Delta_constraint[4][3], 0, "^", label="Second Added Point",
        color="xkcd:teal")

ax.set_ylabel('$F_0$', fontsize=15)
ax.set_xlabel('$\Delta$', fontsize=15)
# ax.set_title('Time for Second Iteration',fontsize=18)

ax.legend()
plt.savefig("outer_plots_2.png",dpi=300)

plt.clf()
fig, ax = plt.subplots()
ax.set_xlim(1.04248,1.04251)
ax.set_ylim(-0.0000000001,0.000000001)

Delta_2 = numpy.arange(1.042,1.043,0.000001)
f1_2 = 1 + numpy.power(Delta_2,4)
f2_2 = numpy.power(Delta_2,2) + numpy.power(Delta_2,4)/12

ax.plot(Delta_2, f1_2 - f2_scale[5] * f2_2, linestyle="-.",
        label="Iteration 6", color="xkcd:magenta")
ax.plot(Delta_2, f1_2 - f2_scale[6] * f2_2, linestyle="dotted",
        label="Iteration 7", color="xkcd:tan")
ax.plot(Delta_2, f1_2 - f2_scale[7] * f2_2, linestyle="dashed",
        label="Iteration 8", color="xkcd:mauve")
ax.plot(Delta_2, numpy.zeros(len(Delta_2)), color="xkcd:black")

ax.plot(Delta_constraint[4][3], 0, "^", label="Second Added Point", color="xkcd:teal")
ax.plot(Delta_constraint[7][3], 0, "D", label="Third Added Point", color="xkcd:forest green")

ax.set_ylabel('$F_0$', fontsize=15)
ax.set_xlabel('$\Delta$', fontsize=15)
# ax.set_title('Time for Second Iteration',fontsize=18)

ax.legend()
plt.savefig("outer_plots_3.png",dpi=300)

# Iteration0
# gap=1.1
# points: 0 [0, 0.1, 1, 1.79769e+308]

# Iteration1
# gap=1.1/1024

# Iteration2
# new points: 0 [0, 0.1, 1, 1.044339780356917847849369292480503379774, 1.797693134862315708145274237317043567981e+308]
# gap=1.1/1024

# Iteration3
# gap=1.1/(1024^2)

# Iteration4
# new points: 0 [0, 0.1, 1, 1.042500153473850784545700476498155026252, 1.044339780356917847849369292480503379774, 1.797693134862315708145274237317043567981e+308]
# gap=1.1/(1024^2)

# Iteration4
# gap=1.1/(1024^3)
