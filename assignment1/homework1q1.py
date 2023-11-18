import matplotlib.pyplot as plt
from engine.truss import TrussLinear, TrussQuadratic

def getLinearSolution(E, L, A, qMin, qMax, F, nElem):
    # Linear Shape Functions Solution
    linear = TrussLinear(E, L, A, nElem)
    linear.formKGlobal()
    linear.formFGlobal(qMin, qMax, 0, F)
    linear.formDirichletNeumann()
    linSol = linear.solve()

    linx = linear.getXField()
    linDisp = linear.getDispField()
    linSigma = linear.E * linear.getStrainField()

    return linear.x_node, linSol, linx, linDisp, linSigma

def getQuadraticSolution(E, L, A, qMin, qMax, F, nElem):
    # Quadratic Shape Function Solution
    quadratic = TrussQuadratic(E, L, A, nElem)
    quadratic.formKGlobal()
    quadratic.formFGlobal(qMin, qMax, 0, F)
    quadratic.formDirichletNeumann()
    quadSol = quadratic.solve()

    quadx = quadratic.getXField()
    quadDisp = quadratic.getDispField()
    quadSigma = quadratic.E * quadratic.getStrainField()

    return quadratic.x_node, quadSol, quadx, quadDisp, quadSigma

def plotwithNodes(ax, label, color, shape, alpha, nodes, sol, solX, solU):
    ax.plot(0, 0, color + shape + '-', label = label, linewidth = 2, alpha = alpha)
    ax.plot(nodes, sol, color + shape , linewidth = 2, alpha = alpha)
    ax.plot(solX, solU, color + '-', linewidth = 2, alpha = alpha)

def plot(ax, label, color, alpha, solX, solU):
    ax.plot(solX, solU, color + '-', label = label, linewidth = 2, alpha = alpha)


if __name__ == '__main__':
    nElem = 3
    L = 1000 # [mm]
    A = 75e3 # [mm^2]
    E = 50e3 # [Mpa]
    F = 100e3 # [N]
    qMin = 200 # [N/mm]
    qMax = 200 # [N/mm]
    
    # Linear Shape Function Solutions
    linSol1x = getLinearSolution(E, L, A, qMin, qMax, F, 1)
    linSol2x = getLinearSolution(E, L, A, qMin, qMax, F, 2)
    linSol4x = getLinearSolution(E, L, A, qMin, qMax, F, 4)
    
    # Quadratic Shape Function Solutions
    quaSol1x = getQuadraticSolution(E, L, A, qMin, qMax, F, 1)
    quaSol2x = getQuadraticSolution(E, L, A, qMin, qMax, F, 2)
    quaSol4x = getQuadraticSolution(E, L, A, qMin, qMax, F, 4)

    # Linear Displacement Plots
    fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=True)
    fig.set_size_inches(12.5, 7.5)
    

    plotwithNodes(ax1, 'Linear 1x', 'r', '^', 0.6, linSol1x[0], linSol1x[1], linSol1x[2], linSol1x[3])
    plotwithNodes(ax1, 'Linear 2x', 'b', 'v', 0.5, linSol2x[0], linSol2x[1], linSol2x[2], linSol2x[3])
    plotwithNodes(ax1, 'Linear 4x', 'g', 'o', 0.4, linSol4x[0], linSol4x[1], linSol4x[2], linSol4x[3])
    
    ax1.set(
        title = 'Displacement Field using Linear Elements',
        ylabel = r'Displacement $\mathbf{[mm]}$',
        xlabel = r'x $\mathbf{[mm]}$',
    )

    ax1.legend(loc='best')
    ax1.grid(True)
    fig.savefig('linear_displacement_vs_x.png', dpi=600)

    # Quadratic Displacement Plots
    fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=True)
    fig.set_size_inches(12.5, 7.5)
    

    plotwithNodes(ax1, 'Quadratic 1x', 'r', '^', 0.6, quaSol1x[0], quaSol1x[1], quaSol1x[2], quaSol1x[3])
    plotwithNodes(ax1, 'Quadratic 2x', 'b', 'v', 0.5, quaSol2x[0], quaSol2x[1], quaSol2x[2], quaSol2x[3])
    plotwithNodes(ax1, 'Quadratic 4x', 'g', 'o', 0.4, quaSol4x[0], quaSol4x[1], quaSol4x[2], quaSol4x[3])
    
    ax1.set(
        title = 'Displacement Field using Quadratic Elements',
        ylabel = r'Displacement $\mathbf{[mm]}$',
        xlabel = r'x $\mathbf{[mm]}$',
    )

    ax1.legend(loc='best')
    ax1.grid(True)
    fig.savefig('quadratic_displacement_vs_x.png', dpi=600)


    # Linear Stress Plots
    fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=True)
    fig.set_size_inches(12.5, 7.5)
    

    plot(ax1, 'Linear 1x', 'r', 0.6, linSol1x[2], linSol1x[4])
    plot(ax1, 'Linear 2x', 'b', 0.5, linSol2x[2], linSol2x[4])
    plot(ax1, 'Linear 4x', 'g', 0.4, linSol4x[2], linSol4x[4])
    
    ax1.set(
        title = 'Stress Field using Linear Elements',
        ylabel = r'$\sigma$ $\mathbf{[MPa]}$',
        xlabel = r'x $\mathbf{[mm]}$',
    )

    ax1.legend(loc='best')
    ax1.grid(True)
    fig.savefig('linear_stress_vs_x.png', dpi=600)


    # Quadratic Stress Plots
    fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=True)
    fig.set_size_inches(12.5, 7.5)
    

    plot(ax1, 'Quadratic 1x', 'r', 0.6, quaSol1x[2], quaSol1x[4])
    plot(ax1, 'Quadratic 2x', 'b', 0.5, quaSol2x[2], quaSol2x[4])
    plot(ax1, 'Quadratic 4x', 'g', 0.4, quaSol4x[2], quaSol4x[4])
    
    ax1.set(
        title = 'Stress Field using Quadratic Elements',
        ylabel = r'$\sigma$ $\mathbf{[MPa]}$',
        xlabel = r'x $\mathbf{[mm]}$',
    )

    ax1.legend(loc='best')
    ax1.grid(True)
    fig.savefig('quadratic_stress_vs_x.png', dpi=600)