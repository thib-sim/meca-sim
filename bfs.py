# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 15:10:51 2022

@author: thibs
"""
class Maillon:

	def __init__(self, valeur, precedent=None, suivant=None):
		self.valeur = valeur
		self.precedent = precedent
		self.suivant = suivant


class File:

	def __init__(self):
		self.longueur = 0
		self.debut = None
		self.fin = None

	def enfiler(self, valeur):
		if self.longueur == 0:
			self.debut = self.fin = Maillon(valeur)
		else:
			self.fin = Maillon(valeur, self.fin)
			self.fin.precedent.suivant = self.fin
		self.longueur += 1


	def defiler(self):
		if self.longueur > 0:
			valeur = self.debut.valeur
			if self.longueur > 1:
				self.debut = self.debut.suivant
				self.debut.precedent = None
			else:
				self.debut = self.fin = None
			self.longueur -= 1
		return valeur


	def estVide(self):
		return self.longueur == 0


def direct_child(aa,bb):
    child=[]
    for cle, valeur in aa.items():
        if valeur in bb:
            child.append(cle)
    return child

def descendance(aa,bb):
    child=bb
    descendant=child
    while child!=[]:
        child=direct_child(aa,child)
        for c in child:
            descendant.append(c)
    return descendant

 

def bfs(G,s):
    P = {s :None}
    D = {}
    a=len(G)
    for u in G:
        D[u]=a+1
    D[s]=0
    Q = File()
    Q.enfiler(s)
    while not(Q.estVide()) :
        u = Q.defiler()
        for v in G[u] :
            if v in P : continue
            P[v]=u
            D[v]=D[u]+1
            Q.enfiler(v)
    return P,D


if __name__ == "__main__":
    G = dict()
    G['a'] = ['b','c']
    G['b'] = ['a','d','e']
    G['c'] = ['a','d']
    G['d'] = ['b','c','e']
    G['e'] = ['b','d','f','g']
    G['f'] = ['e','g']
    G['g'] = ['e','f','h']
    G['h'] = ['g']
    
    
    P,D = bfs(G,'a')
    print(P)
