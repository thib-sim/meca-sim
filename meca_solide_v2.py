# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:51:48 2022

@author: thibs
"""
import numpy as np
import math


class Point:
    def __init__(self,x,y,z,name='',masse=0):
        self.name=name
        self.masse=masse
        self.x=x
        self.y=y
        self.z=z
        
        self.vx=0
        self.vy=0
        self.vz=0 

        self.ax=0
        self.ay=0
        self.az=0 

    def get_position(self):
        return np.matrix([self.x,self.y,self.z])

    def get_vitesse(self):
        return np.matrix([self.vx,self.vy,self.vz])

    def get_accel(self):
        return np.matrix([self.ax,self.ay,self.az])

    def set_position(self,vect):
        self.x=vect[0,0]
        self.y=vect[0,1]
        self.z=vect[0,2]
    
    def set_vitesse(self,vect):
        self.vx=vect[0,0]
        self.vy=vect[0,1]
        self.vz=vect[0,2]     
        
    def set_accel(self,vect):
        self.ax=vect[0,0]
        self.ay=vect[0,1]
        self.az=vect[0,2]
        
    def distance(self,pointA):
        return math.sqrt((self.x-pointA.x)**2+(self.y-pointA.y)**2+(self.z-pointA.z)**2)

    def vecteur(self,pointA):
        return np.matrix([self.x-pointA.x,self.y-pointA.y,self.z-pointA.z])


class repere_2D:
    #par rapport au referentiel_galileleen d'origine (0,0,0) et  d'axe (1,0,0),(0,1,0),(0,0,1)
    
    def __init__(self,point):
        self.origine=point
        x0=np.matrix([1,0,0])
        y0=np.matrix([0,1,0])
        z0=np.matrix([0,0,1])
        axes=np.concatenate((x0,y0,z0))
        a=np.concatenate((axes,point.get_position().transpose()),axis=1)
        b=np.matrix([0,0,0,1])
        self.rep=np.concatenate((a,b))

    def get_axes(self):
        return self.rep[0:3,0:3]
        
    def get_origine(self):
        return self.rep[0:3,3]       

    def set_origine(self,point):
        self.origine=point
        axes=self.get_axes()
        a=np.concatenate((axes,point.get_position().transpose()),axis=1)
        b=np.matrix([0,0,0,1])
        self.rep=np.concatenate((a,b))
    
    def chgt_repere(self,t):
        self.rep=t.dot(self.rep)
        a=self.get_origine()
        self.origine.x=a[0,0]
        self.origine.y=a[1,0]
        self.origine.z=a[2,0]

        
class solid:
    def __init__(self):
        self.centre=Point(0,0,0,'origine solide',0)
        self.masse=0
        self.points={}
        self.repere=repere_2D(Point(0,0,0,'origine solide',0))
        self.inertie=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        self.apply_repere()
                
    def apply_repere(self):
        #calcul les points dans les coordonnées du referentiel d'inertie : points moved
        p=np.concatenate((self.centre.get_position().transpose(),np.matrix([1])))
        pp=self.repere.rep.dot(p)[0:3].transpose()
        self.centre_moved=Point(pp[0,0],pp[0,1],pp[0,2],self.centre.name,self.centre.masse)
        
        self.points_moved={}
        for name in set(self.points):
            p=np.concatenate((self.points[name].get_position().transpose(),np.matrix([1])))
            pp=self.repere.rep.dot(p)[0:3].transpose()
            ppp=Point(pp[0,0],pp[0,1],pp[0,2],self.points[name].name,self.points[name].masse)
            self.points_moved[name]=ppp
            
    def append_point(self,x,y,z,name,m):
        #donner les coordonnées dans le repere de l'objet (par rapport au centre de gravité)
        #m masse en kg du point
        print('ajoute le point ' + name)
        self.points[name]=Point(x,y,z,name,m)
        #actualise les parametres du solide
        self.calcul_centre_gravité()
        self.calcul_matrice_inertie()
        self.apply_repere()
    
    def calcul_centre_gravité(self):
        print('calcul_centre_gravité')
        c_i=np.matrix([0,0,0])
        masse=0
        for name in set(self.points):
            masse+=self.points[name].masse
            c_i=c_i+self.points[name].masse*self.points[name].get_position()
        r=c_i/masse
        self.centre.set_position(r)
        self.masse=masse
        
    def calcul_matrice_inertie(self):
        print('calcul_matrice_inertie')
        m_i=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        masse=0
        #op^(u^op)
        for name in set(self.points):
            masse+=self.points[name].masse
            op=self.points[name].vecteur(self.centre)
            uop=np.cross(self.repere.get_axes(),op)
            m_i=m_i+self.points[name].masse*np.cross(op,uop)
        self.inertie=m_i/masse

    def deplace_origine(self,point):
        #change l'origine du repere equivalent à translation sur le solide
        self.repere.set_origine(point)
        self.apply_repere()
    
    def repositionne_point(self,vect_tr):
        #deplace les points du solide dans son referentiel d'origine
        for name in set(self.points):
            new=self.points[name].get_position()+vect_tr
            self.points[name].set_position(new)
        self.apply_repere()
        
    def deplace_solid(self,t):
        #deplace un solide à partir d'une rotation et d'une translation
        self.repere.chgt_repere(t)
        self.apply_repere()
        

class liaison:
    def __init__(self,solid1_name,solid2_name,pt_ct_1,pt_ct_2,type_liaison):
        self.def_liaison_type={"pivot_z":{"omega":np.matrix([0,0,1]),"translate":np.matrix([0,0,0]),"R":np.matrix([1,1,1]),"M":np.matrix([1,1,0])},
                                "pivot_y":{"omega":np.matrix([0,1,0]),"translate":np.matrix([0,0,0]),"R":np.matrix([1,1,1]),"M":np.matrix([1,0,1])},
                                "pivot_x":{"omega":np.matrix([1,0,0]),"translate":np.matrix([0,0,0]),"R":np.matrix([1,1,1]),"M":np.matrix([0,1,1])},
                                "glissiere_z":{"omega":np.matrix([0,0,0]),"translate":np.matrix([0,0,1]),"R":np.matrix([1,1,0]),"M":np.matrix([1,1,1])},
                                "glissiere_y":{"omega":np.matrix([0,0,0]),"translate":np.matrix([0,1,0]),"R":np.matrix([1,0,1]),"M":np.matrix([1,1,1])},
                                "glissiere_x":{"omega":np.matrix([0,0,0]),"translate":np.matrix([1,0,0]),"R":np.matrix([0,1,1]),"M":np.matrix([1,1,1])}
                                }
        
        self.description={'solid1':solid1_name,
                                     'solid2':solid2_name,
                                     'point_de_contact1':pt_ct_1,
                                     'point_de_contact2':pt_ct_2,
                                     'type_liaison':type_liaison,
                                     "position":{'rot_vect':np.matrix([0.,0.,0.]),'translate':np.matrix([0.,0.,0.])},
                                     "vitesse":{'rot_vect':np.matrix([0.,0.,0.]),'translate':np.matrix([0.,0.,0.])}
                                     }
        
        
    def get_type(self):       
        l_t=self.description['type_liaison']
        return self.def_liaison_type[l_t]['omega'],self.def_liaison_type[l_t]['translate']        

    def get_type2(self):       
        l_t=self.description['type_liaison']
        return self.def_liaison_type[l_t]['R'],self.def_liaison_type[l_t]['M']    
        
    def get_solid(self):        
        return self.description['solid1'],self.description['solid2']

    def get_vitesse(self):        
        return self.description['vitesse']['rot_vect'],self.description['vitesse']["translate"]
        
    def set_vitesse(self,norm_omega,norm_translate):    
        omega,translate=self.get_type()
        self.description['vitesse']['rot_vect']=omega*norm_omega
        self.description['vitesse']["translate"]=norm_translate*translate
        
    def get_position(self):        
        return self.description['position']['rot_vect'],self.description['position']["translate"]
        
    def set_position(self,norm_omega,norm_translate):    
        type_liaison=self.description["type_liaison"]
        omega=norm_omega*self.def_liaison_type[type_liaison]['omega']
        translate=norm_translate*self.def_liaison_type[type_liaison]['translate']
        self.description['position']['rot_vect']=omega
        self.description['position']["translate"]=translate
        
    def update_position(self,omega,translate):    
        self.description['position']['rot_vect']+=omega
        self.description['position']["translate"]+=translate
        
    def update_vitesse(self,omega,translate):    
        self.description['vitesse']['rot_vect']+=omega
        self.description['vitesse']["translate"]+=translate
        
    def increase_vit(self):
        norm_omega=2*math.pi/60
        norm_translate=10
        omega,translate=self.get_type()
        self.description['vitesse']['rot_vect']+=omega*norm_omega
        self.description['vitesse']["translate"]+=norm_translate*translate

    def decrease_vit(self):
        norm_omega=2*math.pi/60
        norm_translate=10
        omega,translate=self.get_type()
        self.description['vitesse']['rot_vect']-=omega*norm_omega
        self.description['vitesse']["translate"]-=norm_translate*translate


class cinematique:
    def __init__(self):
        self.t_cinematiques={}
        self.t_result=torseur_cinematique()
        
    def append_torseur(self,name,solid,pt_application,repere_name,omega,translate):
        print('ajoute torseur ' + name)
        # repere = soit ref d' inertie soit piece
        # pt d'application doit etre exprimé dans le referentiel de la piece ex : pour le centre de gravité choisir solid.centre
        t_d2=torseur_cinematique()
        t_d2.init_from_solid(solid,pt_application,repere_name)
        t_d2.set_point(pt_application,translate)
        t_d2.set_rslt(omega)
        self.t_cinematiques[name]=t_d2
        self.compute_resultante()
    
    def compute_resultante(self):
        #sous entends une chaine cinematique ouverte avec la premiere liaison relié au bati
        #calcul le mvt du end effecteur par rapport au ref d'inertie
        print("calcul le torseur vitesse")
        point_ref=self.t_result.point_ref
        self.t_result.set_rslt()
        self.t_result.set_moment()
        for nom in self.t_cinematiques:
            td=self.t_cinematiques[nom]
            if td.repere_name=="inertie":
                td.chge_point(point_ref)
            elif td.repere_name=="piece":
                td.chge_repere(td.solid)
                td.chge_point(point_ref)
            else:
                print("WARNING : repere_name non conforme")
            self.t_result=self.t_result.compose(td)
            
    def list_t(self):
        return set(self.t_cinematiques)

    def get_resultante(self):
        return self.t_result.get_rslt()
    
    def get_moment(self):
        return self.t_result.get_moment()
  
    def change_resultante_point(self,point):
        print('update calcul force')
        # repere = soit ref d' inertie soit piece
        # pt d'application doit etre exprimé dans le referentiel de la piece ex : pour le centre de gravité choisir solid.centre
        self.t_result.chge_point(point)
        self.compute_resultante()  
        
    def get_norm_vit(self):
        return math.sqrt(self.get_moment().dot(self.get_moment().transpose()))
        
class force:
    def __init__(self):
        #la force s'applique :
        #    sur un point du solide si la force est ponctuel
        #    sur une face du solide si la force est surfacique
        #    sur le centre de gravité si la force est volumique
        self.forces={}
        self.t_d_result=torseur_dynamique()
    
    def list_force(self):
        return set(self.forces)

    def get_resultante(self):
        return self.t_d_result.get_rslt()
    
    def get_moment(self):
        return self.t_d_result.get_moment()
    
    def append_force(self,name,vect,solid,pt_application,repere_name):
        print('ajoute la force ' + name)
        # repere = soit ref d' inertie soit piece
        # pt d'application doit etre exprimé dans le referentiel de la piece ex : pour le centre de gravité choisir solid.centre
        t_d2=torseur_dynamique()
        t_d2.init_from_solid(solid,pt_application,repere_name)
        t_d2.set_rslt(vect)
        self.forces[name]=t_d2
        self.compute_resultante()    
        
    def compute_resultante(self):
        print("calcul_force_resultante dans le ref d'inertie")
        point_ref=self.t_d_result.point_ref
        self.t_d_result.set_rslt()
        self.t_d_result.set_moment()
        for nom in self.forces:
            td=self.forces[nom]
            if td.repere_name=="inertie":
                td.chge_point(point_ref)
            elif td.repere_name=="piece":
                td.chge_repere(td.solid)
                td.chge_point(point_ref)
            else:
                print("WARNING : repere_name non conforme")
            self.t_d_result=self.t_d_result.compose(td)


    def change_resultante_point(self,point):
        print('update calcul force')
        # repere = soit ref d' inertie soit piece
        # pt d'application doit etre exprimé dans le referentiel de la piece ex : pour le centre de gravité choisir solid.centre
        self.t_d_result.chge_point(point)
        self.compute_resultante()    


class transform:
    #cette class genere une matrice de transformation homogene à partir des angles d'euler et d'un vecteur de translation
    def __init__(self,rot_v=np.matrix([0,0,0]),t_vect=np.matrix([0,0,0]),theta=0.):
        rot_mat=rotation_matrix()
        rot_mat.create_from_vect(rot_v,theta)
        self.rotation=rot_mat.r
        self.translate=t_vect.transpose()
        self.create_transform()
        self.create_inverse_transform()
        
    def create_transform(self):
        a=np.concatenate((self.rotation,self.translate),axis=1)
        b=np.matrix([0,0,0,1])
        self.transform=np.concatenate((a,b))
        
    def get_rotation(self):
        return self.transform[0:3,0:3]

    def get_translate(self):
        return self.transform[0:3,3]
        
    def create_inverse_transform(self):
        print('')
        a=self.get_rotation()
        b=self.get_translate()
        c=a.transpose()
        d=-c.dot(b)
        e=np.concatenate((c,d),axis=1)
        f=np.matrix([0,0,0,1])
        self.inverse_transform=np.concatenate((e,f))        

    def get_inv_rotation(self):
        return self.inverse_transform[0:3,0:3]

    def get_inv_translate(self):
        return self.inverse_transform[0:3,3]        
    
    def input_rot_mat(self,rot_mat):
        a=np.concatenate((rot_mat,self.translate),axis=1)
        b=np.matrix([0,0,0,1])
        self.transform=np.concatenate((a,b))
    
        
class rotation_matrix:
    #cette class genere une matrice de rotation à partir des 3 angles d'euler
    def __init__(self):
        self.r=np.eye(3)
        
    def create_from_vect(self,vect,theta):
        u_x=vect[0,0]
        u_y=vect[0,1]
        u_z=vect[0,2]    

        c=math.cos(theta)
        s=math.sin(theta)
        
        a11=u_x*u_x*(1-c)+c
        a12=u_x*u_y*(1-c)-u_z*s
        a13=u_x*u_z*(1-c)+u_y*s
        a21=u_x*u_y*(1-c)+u_z*s
        a22=u_y*u_y*(1-c)+c
        a23=u_y*u_z*(1-c)-u_x*s
        a31=u_x*u_z*(1-c)-u_y*s
        a32=u_y*u_z*(1-c)+u_x*s
        a33=u_z*u_z*(1-c)+c

        mat=np.matrix([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
        self.r=mat
    
    def create_from_euler(self,rot_vect=np.matrix([0,0,0])):
        angle_x=rot_vect[0,0]
        angle_y=rot_vect[0,1]
        angle_z=rot_vect[0,2]
        
        rotation_z = np.matrix([
            [math.cos(angle_z), -math.sin(angle_z), 0],
            [math.sin(angle_z), math.cos(angle_z), 0],
            [0, 0, 1],
        ])
        
        rotation_y = np.matrix([
            [math.cos(angle_y), 0, math.sin(angle_y)],
            [0, 1, 0],
            [-math.sin(angle_y), 0, math.cos(angle_y)],
        ])
        
        rotation_x = np.matrix([
            [1, 0, 0],
            [0, math.cos(angle_x), -math.sin(angle_x)],
            [0, math.sin(angle_x), math.cos(angle_x)],
        ])
        
        c=np.dot(rotation_y,rotation_z)
        rotation_matrix=np.dot(rotation_x,c)
        self.r=rotation_matrix
        


class torseur:
    def __init__(self):
        #repere_name = inertie
        #"creation d'un torseur dynamique")
        #"moment dynamique = au point fixe")
        self.torseur=np.matrix([[0.,0.,0.],[0.,0.,0.]])
        self.repere_name='inertie'
        self.point_ref=Point(0,0,0)
        self.solid=0
        
    def init_from_solid(self,solid,pt,repere_name):
        #pt  : nom d'un pt du solid exprimé dans le repere choisi => choisir pt_moved si inertie
        #solid : objet solide appliquant la force
        #repere_name : soit inertie soit piece
        self.repere_name=repere_name
        self.solid=solid
        self.point_ref=pt
        
    def set_rslt(self,p=np.matrix([0,0,0])):
        self.torseur[0,0:3]=p
        
    def get_rslt(self):
        return self.torseur[0,0:3]     

    def get_moment(self):
        return self.torseur[1,0:3]     
    
    def set_moment(self,p=np.matrix([0,0,0])):
        self.torseur[1,0:3]=p

    def set_point(self,pointA,mA=np.matrix([0,0,0])):
        #faire attention que le point A soit defini ds le meme repere que le torseur
        self.torseur[1,0:3]=mA    
        self.point_ref=pointA

    def chge_point(self,pointB):
        #faire attention que le point B soit defini ds le meme repere que celui qu il remplace
        pointA=self.point_ref
        vect=(pointB.get_position()-pointA.get_position())
        vB=self.get_moment()+np.cross(-vect,self.get_rslt())
        self.set_point(pointB,vB)

    def chge_repere(self,solid):
        self.solid=solid
        #passer les composantes du torseur ds le nouveau repere
        #exprimer le point dans le nouveau repere
        if self.repere_name=='inertie':
            
            self.repere_name='piece'
            axes=self.solid.repere.get_axes()
            axes_t=axes.transpose()
            abs_rslt=axes_t.dot(self.get_rslt().transpose()).transpose()
            self.set_rslt(abs_rslt)
            abs_monent=axes_t.dot(self.get_moment().transpose()).transpose()
            
            t=transform()
            t.transform=self.solid.repere.rep
            t.create_inverse_transform()
            t_inv=t.inverse_transform
            p=np.concatenate((self.point_ref.get_position().transpose(),np.matrix([1])))
            pp=t_inv.dot(p)[0:3].transpose()
            point_ref2=Point(pp[0,0],pp[0,1],pp[0,2],self.point_ref.name,self.point_ref.masse)
            
            self.set_point(point_ref2,abs_monent)
            
        elif self.repere_name=='piece':
            
            self.repere_name='inertie'  
            axes=self.solid.repere.get_axes()
            abs_rslt=axes.dot(self.get_rslt().transpose()).transpose()
            self.set_rslt(abs_rslt)
            abs_monent=axes.dot(self.get_moment().transpose()).transpose()
                        
            p=np.concatenate((self.point_ref.get_position().transpose(),np.matrix([1])))
            pp=self.solid.repere.rep.dot(p)[0:3].transpose()
            point_ref2=Point(pp[0,0],pp[0,1],pp[0,2],self.point_ref.name,self.point_ref.masse)
            
            self.set_point(point_ref2,abs_monent)
            
        else:
            print('unknown referentiel')

    def compose(self,tenseur):
        #il faut etre dans le meme referentiel au même point
        print('')
        if self.repere_name!='inertie' or tenseur.repere_name!='inertie':
            print("WARNING :les tenseurs ne sont pas exprimés dans le repere d'inertie")
        else:
            if self.point_ref.name!=tenseur.point_ref.name:
                print('WARNING :les tenseurs ne sont pas exprimés au meme point')
            else:
                a=torseur()
                a.torseur=self.torseur+tenseur.torseur
                a.point_ref=self.point_ref
                return a


class torseur_dynamique(torseur):
    def __init__(self):
        super().__init__()


class torseur_cinematique(torseur):
    def __init__(self):
        super().__init__()


class torseur_cinetique(torseur):
    def __init__(self):
        super().__init__()
        
      
if __name__ == "__main__":
    
    a=solid()
    a.append_point(0,50,0,'rar_g',10)
    a.append_point(0,-50,0,'rar_d',10)
    a.append_point(80,50,0,'rav_g',10)
    a.append_point(80,-50,0,'rav_d',10)
    a.apply_repere()
    print(a.inertie)
    
    t=transform(np.matrix([0,0,1]),t_vect=np.matrix([0,0,0]),theta=math.pi/2)
    a.deplace_solid(t.transform)    
    
    f=force()
    f.append_force('mot_g', np.matrix([10,0,0]), a, a.points['rar_g'],'piece')
    f.append_force('mot_d', np.matrix([10,0,0]), a,  a.points['rar_d'],'piece')    
    f.compute_resultante()
    r=rotation_matrix()

    t_c=torseur_cinematique()
    t_c.set_rslt(np.matrix([0,0,10]))
    t_c.get_rslt()
    t_c.get_moment()
    t_c.chge_point(Point(20,0,0))
    t_c.get_moment()

    t_d=torseur_dynamique()
    t_d.init_from_solid(a, a.points['rar_g'],'piece')
    t_d.set_rslt(np.matrix([10,0,0]))
    t_d.get_rslt()
    t_d.chge_point(a.points['rar_d'])
    t_d.chge_repere(t_d.solid)

    t_d2=torseur_dynamique()
    t_d2.init_from_solid(a, a.points['rar_d'],'piece')
    t_d2.set_rslt(np.matrix([-10,0,0]))
    t_d2.get_rslt()
    t_d2.chge_repere(t_d.solid)
       
    b=t_d.compose(t_d2)
    np.linalg.inv(a.inertie)
