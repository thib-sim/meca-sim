# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 12:52:31 2022

@author: thibs
"""
from meca_solide_v2 import *
from threading import Thread
import time
from bfs import bfs,descendance

class system_meca:
    def __init__(self):
        self.run=False
        self.temps=0
        self.dt=0.1
        self.solids={}
        self.liaisons={}
        self.chaine_cinematique={}
        self.force_ext=force()
        self.repere=repere_2D(Point(0,0,0,'origine solide',0))
        self.torseurs_cinetique={}
        self.moteurs={}
        self.cinematique=cinematique()

    def chaine(self):
        chaine={}
        liaison=[]
        for l in self.liaisons:
            solid1=self.liaisons[l].description['solid1']
            solid2=self.liaisons[l].description['solid2']
            if chaine=={}:
                chaine[solid1]=[solid2]
                chaine[solid2]=[solid1]
                
            else:
                if solid1 in chaine:
                    if solid2 not in chaine[solid1]:
                        chaine[solid1].append(solid2)
                else:
                    chaine[solid1]=[solid2]
                    
                if solid2 in chaine:
                    if solid1 not in chaine[solid2]:
                        chaine[solid2].append(solid1)
                else:
                    chaine[solid2]=[solid1]
                    
        P,D=bfs(chaine,'bati')
        return P

    def create_moteurs(self):
        #create torseur with unitary componant
        for liaison_name in self.liaisons:
            liaison=self.liaisons[liaison_name]
            pt_ct_name=liaison.description['point_de_contact2']
            solid=liaison.description['solid2']
            pt=self.solids[solid].points_moved[pt_ct_name]
            axes=self.solids[solid].repere.get_axes()
            #utilisation de la forme des torseurs cinematiques
            M,R=liaison.get_type()
            abs_R=axes.dot(R.transpose()).transpose()
            abs_M=axes.dot(M.transpose()).transpose()
            a=torseur_dynamique()
            a.set_point(pt.get_position(),abs_M)
            a.set_rslt(abs_R)
            self.moteurs[liaison_name]=a
        
    def create_torseur_cinetique(self):
        #create torseur with unitary componant
        for liaison_name in self.liaisons:
            liaison=self.liaisons[liaison_name]
            pt_ct_name=liaison.description['point_de_contact2']
            solid=liaison.description['solid2']
            pt=self.solids[solid].points_moved[pt_ct_name]
            axes=self.solids[solid].repere.get_axes()
            R,M=liaison.get_type2()
            abs_R=axes.dot(R.transpose()).transpose()
            abs_M=axes.dot(M.transpose()).transpose()
            a=torseur_cinetique()
            a.set_point(pt.get_position(),abs_M)
            a.set_rslt(abs_R)
            self.torseurs_cinetique[liaison_name]=a

    def calcul_torseur_cinematique(self):
        #create torseur with unitary componant
        self.cinematique=cinematique()
        for liaison_name in self.liaisons:
            liaison=self.liaisons[liaison_name]
            pt_ct_name=liaison.description['point_de_contact2']
            solid_name=liaison.description['solid2']
            solid=self.solids[solid_name]
            pt=solid.points[pt_ct_name]
            omega,translate=liaison.get_vitesse()
            self.cinematique.append_torseur(solid_name,solid,pt,'piece',omega,translate)
       
        
    def apply_gravity(self):
        for  solid_name in self.solids:
            print(solid_name)
            solid=self.solids[solid_name]
            self.force_ext.append_force('poids '+solid_name, np.matrix([0,0,-10*solid.masse]), solid, solid.centre_moved, 'inertie')

    def apply_bati(self,solid_name):
            solid=self.solids[solid_name]
            x=self.force_ext.get_resultante()
            print(x)
            self.force_ext.append_force('resultante bati', -x, solid, solid.centre_moved, 'inertie')

            
    def append_solid(self,solid,solid_name):
        #ajoute un solid
        print('ajoute le solid ' + solid_name)
        self.solids[solid_name]=solid

    def append_liaison(self,liaison_name,solid1_name,solid2_name,pt_ct_1,pt_ct_2,type_liaison):
        #donner les noms des points de contacts de chacun des solides
        print('ajoute la liaison entre ' + solid1_name+' and '+ solid2_name)
        #repositionne le repere du solid 2 sur le pt de contact du solid 2
        #puis positionne ce point sur point de contact du solide1
        old_origine=self.solids[solid2_name].repere.origine
        new_origine=self.solids[solid2_name].points[pt_ct_2]
        vect=old_origine.vecteur(new_origine)
        self.solids[solid2_name].repositionne_point(vect)
        coord=self.solids[solid1_name].points_moved[pt_ct_1]
        self.solids[solid2_name].deplace_origine(coord) 
        ll=liaison(solid1_name,solid2_name,pt_ct_1,pt_ct_2,type_liaison)
        self.liaisons[liaison_name]=ll
        self.chaine_cinematique=self.chaine()
        
    def set_vitesse_liaison(self,liaison_name,norm_omega,norm_translate):
        self.liaisons[liaison_name].set_vitesse(norm_omega,norm_translate)
 
    def compute_mvt_to_apply_from_vit(self):
        dt=self.dt
        self.list_mvt=[]
        for liaison_name in self.liaisons:
            rot_vect,translate= self.liaisons[liaison_name].get_vitesse()
            type_omega,type_translate=self.liaisons[liaison_name].get_type()
            norm_omega=(rot_vect*type_omega.transpose())*dt
            norm_translate=(translate*type_translate.transpose())*dt
            
            self.list_mvt.append({'liaison_name':liaison_name,'norm_omega':norm_omega,'norm_translate':norm_translate})
        
        
    def move_system(self,liaison_name,norm_omega,norm_translate):
        #calcul la transformation puis l'applique au systeme
    
        solid1,solid2=self.liaisons[liaison_name].get_solid()
  
        vect=self.solids[solid2].repere.get_origine().transpose()        
        axes=self.solids[solid1].repere.get_axes()
        
        omega,translate0=self.liaisons[liaison_name].get_type()

        translate=norm_translate*translate0
        self.liaisons[liaison_name].update_position(norm_omega*omega,translate)
        
        #calcul du vecteur translation à parir de l'orientation de la piece
        abs_translate=axes.dot(translate.transpose()).transpose()
        #calcul l'axe de rotation à parir de l'orientation de la piece
        abs_omega=axes.dot(omega.transpose()).transpose()
        #calcul de la  transform de la piece 2 par rapport  à la  piece 1        
        t1=transform(abs_omega,abs_translate,norm_omega)
        tt=t1.transform
        
        #chgt de repere de la transform : expression par rapport au  ref d'inertie
        #ppp=p_inv*tt*p avec matrice de chgt de repere qui permet d'appliquer la transorm dans le ref d'inertie (0,x,y,z)=ref du bati
        p=transform(np.matrix([0,0,0]),-vect).transform
        p_inv=transform(np.matrix([0,0,0]),vect).transform
        pp=tt.dot(p)
        ppp=p_inv.dot(pp)

        #parcours de la chaine cinematique et applique le mvt

        bb=[solid2]
        aa=self.chaine_cinematique
        descendant=descendance(aa,bb)
        
        for i in descendant:
            self.solids[i].deplace_solid(ppp)

          
    def apply_list_mvt(self):
        for mvt in self.list_mvt:
            self.move_system(mvt['liaison_name'],mvt['norm_omega'],mvt['norm_translate'])
        self.list_mvt=[]
        self.temps+=self.dt
        
    def run_simu(self):
        self.run=True
        while self.run:
            self.compute_mvt_to_apply_from_vit()
            self.apply_list_mvt()
            self.apply_gravity()
            self.apply_bati('bati')
            self.calcul_torseur_cinematique()
            print('force')
            print(self.force_ext.get_resultante())
            print(self.force_ext.get_moment())
            
            print('point_ref')
            self.cinematique.change_resultante_point(self.solids['bras1'].points_moved['a'])
            print(self.cinematique.t_result.point_ref.get_position())
            print('vitesse')
            print(self.cinematique.get_resultante())
            print(self.cinematique.get_moment())
            print(self.cinematique.get_norm_vit())
            print(self.temps)
            time.sleep(0.1)
    
    def acquisition_commute(self):
        self.run= not self.run
        if self.run:
            self.acq_threads=Thread(target=self.run_simu)
            self.acq_threads.start()
        else:
            self.acq_threads.join()
            print('fin de la simulation')    
    


if __name__ == "__main__":
    a=solid()
    a.append_point(400,400,0,'g',10)    
    a.append_point(400,400,0,'h',10)
    
    b=solid()
    b.append_point(10,100,0,'a',10)
    b.append_point(10,50,0,'b',10)
    b.append_point(10,80,0,'a_bis',10)
    
    c=solid()
    c.append_point(60,0,0,'c',10)
    c.append_point(10,0,0,'d',10)


    d=solid()
    d.append_point(0,0,0,'e',10)
    d.append_point(0,50,0,'f',10)

    e=solid()
    e.append_point(0,0,0,'i',10)
    e.append_point(0,50,0,'j',10)
    
    
    thib=system_meca()
    thib.append_solid(a,'bati')
    thib.append_solid(b,'bras1')
    thib.append_solid(c,'bras2')
    thib.append_solid(d,'bras3')
    thib.append_solid(e,'bras4')
    thib.append_liaison('pivot01', 'bati', 'bras1', 'h', 'b','pivot_z')
    thib.append_liaison('pivot12', 'bras1', 'bras2', 'a', 'd','pivot_z')
    thib.append_liaison('pivot23', 'bras2', 'bras3', 'c', 'e','glissiere_y')
    thib.append_liaison('pivot14', 'bras1', 'bras4', 'a_bis', 'i','pivot_z')
    print(thib.chaine_cinematique)
    
    thib.set_vitesse_liaison('pivot12',2*math.pi/60,0)
    thib.compute_mvt_to_apply_from_vit()
    thib.apply_list_mvt()
    thib.apply_gravity()
    thib.apply_bati('bras1')
    thib.create_torseur_cinetique()
    thib.calcul_torseur_cinematique()
    thib.create_moteurs()

    print(thib.force_ext.get_moment())
    
    thib.set_vitesse_liaison('pivot12',2*math.pi/60,0)
    thib.compute_mvt_to_apply_from_vit()
    thib.apply_list_mvt()    
    
    thib.apply_gravity()
    print(thib.force_ext.get_moment())


    aa=thib.chaine_cinematique
    
    bb=['bras1']
    

    
