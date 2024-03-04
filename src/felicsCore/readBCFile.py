def initBCdict(param, lineIDs, SolutionDict):

	bcDict = {}
	for k in param.SolutionList:
		if k == 'u':
			for kk in param.VelocityComponents:
				bcDict[k+kk] = []
				for kkk in lineIDs:
					bcDict[k+kk].append({'ID':kkk, 'type': 'Dirichlet', 'value':0})
		else:
			bcDict[k] = []
			for kkk in lineIDs:
				bcDict[k].append({'ID':kkk, 'type': 'Dirichlet', 'value':0})

	return bcDict

def readBCFile(param, FEMSpaces):
	from functions import printDebug
	from fenics import DirichletBC,MeshFunction
	file = open(param.BCs.BCsFilePath,'r')
	bcDict = eval(file.read())
	bcDict = bcDict
	file.close()
	BClist = []
	mesh=FEMSpaces.P2.mesh()
	boundaries = MeshFunction('size_t', mesh, param.Case.MeshFilePath.split('.')[0] + '_facet_region.xml')
	for k,m in zip(list(bcDict.keys()),range(0,len(bcDict.keys()))):
		# Get index of equation/variable i_eqn and if needed the index of the velocity component
		if k[0]=='u' and k[1] in param.VelocityComponents:
			i_eqn=0
			i_component=param.VelocityComponents.index(k[1])
		elif k in param.SolutionList:
			i_eqn = param.SolutionList.index(k)
		else:
			continue

		for mm in range(0,len(bcDict[k])) :
			if bcDict[k][mm]['type']=='Dirichlet':


				printDebug(param.debug,"Adding Dirichlet BC for "+str(k)+ " in equation "+str(i_eqn)+" with value "+str(bcDict[k][mm]['value'])+" on boundary with index "+str(bcDict[k][mm]['ID']))
				if k in ['u'+ component for component in param.VelocityComponents]:
					BClist.append( DirichletBC(FEMSpaces.VMixed.sub(i_eqn).sub(i_component), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )


				else:
					BClist.append( DirichletBC(FEMSpaces.VMixed.sub(i_eqn), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )
	return BClist, bcDict
