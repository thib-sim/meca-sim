# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:45:39 2022

@author: thibs
"""
from system_meca import *
import tkinter as tk



class affichage_system_meca:
    def __init__(self,win,system_meca):
        self.win=win
        self.height=800
        self.width=800
        self.win.grid()
        self.win.protocol("WM_DELETE_WINDOW", self.quit_app)
        self.system_meca=system_meca
        self.cadre_menu=tk.Frame(self.win)
        self.cadre_menu.grid(row=1,sticky='NW')
        self.cadre_graph=tk.Frame(self.win)
        self.cadre_graph.grid(row=2)
        self.dessine()
        self.affichage_liaison()
        self.system_meca.acquisition_commute()
        self.suite_mvt()
        
    def dessine(self):
        #print('dessine')
        self.canvas=tk.Canvas(self.cadre_graph,height=self.height,width=self.width)
        self.canvas.grid(row=0,column=0)
        self.canvas.delete('ALL')
        h=6

        for solid_name in self.system_meca.solids:
            solid=self.system_meca.solids[solid_name]
            points=solid.points_moved
            centre=solid.centre_moved
            for name in set(points):
                self.canvas.create_rectangle(int(points[name].x)-int(h/2), int(points[name].y)-int(h/2), int(points[name].x)+int(h/2), int(points[name].y)+int(h/2), fill='red', outline='')
                self.canvas.create_line(int(points[name].x), int(points[name].y),int(centre.x), int(centre.y), width=1)
    
    def affichage_liaison(self):
        i=0
        self.graph_button={}
        for liaison_name in self.system_meca.liaisons:
            self.graph_button[liaison_name]=tk.Button(self.cadre_menu,text="vit +",command=self.system_meca.liaisons[liaison_name].increase_vit)
            self.graph_button[liaison_name].grid(row=0+i,column=1)
            self.graph_button[liaison_name]=tk.Button(self.cadre_menu,text="vit -",command=self.system_meca.liaisons[liaison_name].decrease_vit)
            self.graph_button[liaison_name].grid(row=0+i,column=2)
            i+=1

    

    def suite_mvt(self):
        self.dessine()
        self.win.after(200,self.suite_mvt)
        
        
    def quit_app(self):
        print('cloture acquisition')
        self.system_meca.acquisition_commute()
        self.win.destroy()

if __name__ == "__main__":
    
    a=solid()
    a.append_point(400,400,0,'g',10)    
    a.append_point(400,400,0,'h',10)
    
    b=solid()
    b.append_point(0,100,0,'a',10)
    b.append_point(0,0,0,'b',10)

    
    c=solid()
    c.append_point(50,0,0,'c',10)
    c.append_point(0,0,0,'d',10)


    d=solid()
    d.append_point(0,0,0,'e',10)
    d.append_point(0,50,0,'f',10)

    e=solid()
    e.append_point(0,0,0,'g',10)
    e.append_point(0,50,0,'h',10)

    f=solid()
    f.append_point(0,0,0,'j',10)
    f.append_point(0,50,0,'k',10)


    thib=system_meca()
    thib.append_solid(a,'bati')
    thib.append_solid(b,'bras1')
    
    thib.append_solid(c,'bras2')
    thib.append_solid(d,'bras3')
    thib.append_solid(e,'brasaux')
    thib.append_solid(f,'brasaux2')
    
    thib.append_liaison('pivot01', 'bati', 'bras1', 'h', 'b','pivot_z')
    thib.append_liaison('pivot1aux', 'bras1', 'brasaux', 'a', 'g','glissiere_y')    
    thib.append_liaison('pivot12', 'bras1', 'bras2', 'a', 'd','pivot_y')
    thib.append_liaison('pivot23', 'bras2', 'bras3', 'c', 'e','pivot_x')
    thib.append_liaison('pivotauxaux2', 'brasaux', 'brasaux2', 'h', 'j','pivot_x')
    
 
    thib.set_vitesse_liaison('pivot12',2*math.pi/60,0)
    
    thib.set_vitesse_liaison('pivot01',2*math.pi/60,0)
    thib.force_ext.change_resultante_point(a.centre_moved)
    thib.cinematique.change_resultante_point(b.points_moved['a'])
    win=tk.Tk()
    ap=affichage_system_meca(win,thib)
    win.mainloop()
    
    
    """    
    thib.move_system('pivot12',-math.pi/8,0)
    print(thib.solids["bras3"].points_moved["f"].get_position())
    win=tk.Tk()
    ap=affichage_system_meca(win,thib)
    win.mainloop()
    
    thib.move_system('pivot23',-math.pi/2,100)
    print(thib.solids["bras3"].points_moved["f"].get_position())
    win=tk.Tk()
    ap=affichage_system_meca(win,thib)
    win.mainloop()

    thib.move_system('pivot12',math.pi/2,0)
    print(thib.solids["bras3"].points_moved["f"].get_position())
    
    win=tk.Tk()
    ap=affichage_system_meca(win,thib)
    win.mainloop()

    thib.move_system('pivot23',math.pi/2,-100)
    print(thib.solids["bras3"].points_moved["f"].get_position())
    
    win=tk.Tk()
    ap=affichage_system_meca(win,thib)
    win.mainloop()
    """

