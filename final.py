import math
import matplotlib.pyplot as plt
import numpy as np
import jinja2
import weasyprint

class Section:
    def __init__(self, section_type, E, Fy, bf_top, bf_bot, tf_top, tf_bot, hw, tw):
        self.section_type = section_type
        self.bf_top = bf_top
        self.bf_bot = bf_bot
        self.tf_top = tf_top
        self.tf_bot = tf_bot
        self.hw = hw
        self.tw = tw

        self.h = hw+tf_bot+tf_top

        self.Ag = bf_top*tf_top + bf_bot*tf_bot + hw*tw

        self.yg_from_top = ((bf_bot*(tf_top**2)*0.5)+(bf_bot*tf_bot*(self.h-0.5*tf_bot))+(hw*tw*(tf_top+0.5*hw)))/self.Ag
        
        self.I_top = (bf_top*tf_top**3)/12 + bf_top*tf_top*(self.yg_from_top-0.5*tf_top)**2 + ((self.yg_from_top-tf_top)**3*tw)/12 + (self.yg_from_top-tf_top)*tw*(self.yg_from_top-0.5*(self.yg_from_top+tf_top))**2
        self.I_bot = (bf_bot*tf_bot**3)/12 + bf_bot*tf_bot*(self.h-0.5*tf_bot-self.yg_from_top)**2 + ((self.h-tf_bot-self.yg_from_top)**3*tw)/12 + (self.h-tf_bot-self.yg_from_top)*tw*(0.5*(self.yg_from_top+self.h-tf_bot)-self.yg_from_top)**2
        self.I = self.I_top+self.I_bot
        
        self.Sx_top = self.I/self.yg_from_top
        self.Sx_bot = self.I/(self.h-self.yg_from_top)
        
        x = (bf_bot*tf_bot - bf_top*tf_top + hw*tw)/(2*tw)
        if x<=hw:
            if x>=0:
                self.Z = bf_top*tf_top*(x+0.5*tf_top) + x*tw*0.5*x + tw*0.5*(hw-x)**2 + bf_bot*tf_bot*(hw-x+0.5*tf_bot)
            else:
                x = (hw*tw + bf_bot*tf_bot +tf_top*bf_top)/(2*bf_top)
                c = tf_top - x
                self.Z = bf_top*0.5*x**2 + bf_top*0.5*c**2 + hw*tw*(c+0.5*hw) + bf_bot*tf_bot*(c+hw+0.5*tf_bot)
        else:
            x = (hw*tw + bf_top*tf_top +tf_bot*bf_bot)/(2*bf_bot)
            c = tf_bot - x
            self.Z = bf_bot*0.5*x**2 + bf_bot*0.5*c**2 + hw*tw*(c+0.5*hw) + bf_top*tf_top*(c+hw+0.5*tf_top)

        self.landa = self.hw/self.tw
    
        self.landa_p_web = 3.76 * math.sqrt(E/Fy)
        self.landa_r_web = 5.7 * math.sqrt(E/Fy)

        

    def print(self):
        text = "h = {:0.2f}\n<br>Ag = {:0.2f}\n<br>yg_from_top = {:0.2f}\n<br>I_top = {:0.2f}\n<br>I-bot = {:0.2f}\n<br>I = {:0.2f}<br>Sx_top = {:0.2f}\n<br>Sx_bot = {:0.2f}\n<br>Z = {:0.2f}\n<br>landa_p_web = {:0.2f}\n<br>landa_r_web = {:0.2f}\n<br>landa_web = {:0.2f}".format(self.h,self.Ag,self.yg_from_top,self.I_top,self.I_bot,self.I,self.Sx_top,self.Sx_bot,self.Z,self.landa_p_web,self.landa_r_web,self.landa)
        return text

class SFunction:
    def __init__(self,start,factor,power):
        self.start= start
        self.factor= factor
        self.power= power
    
    def valueAt(self, x):
        if x <= self.start:
            return 0
        else:
            return self.factor * ((x-self.start)**self.power)
    
    def integrate(self,number_of_integrations = 1):
        factor = self.factor
        power = self.power
        for i in range(number_of_integrations):
            factor = factor/(power + 1)
            power = power + 1
        
        return SFunction(self.start, factor, power)
    
    def print(self):
        t = "{:0.2f}*[x-{:0.2f}]^{:0.2f} ".format(self.factor,self.start,self.power)
        return t

class Reaction:
    def __init__(self,cord):
        self.cord = cord
        self.sf_in_shear = SFunction(cord,1,0)
        self.sf_in_deflection = self.sf_in_shear.integrate(3)
        self.sf_in_moment = self.sf_in_shear.integrate()
        self.value = 0
    
    def print(self):
        t = "x = {:0.2f} , value = {:0.2f}<br>".format(self.cord,self.value)
        return t

def sectioncalc(E,Fy):
    section_type = str(input("Enter the section type I or U: "))
    if section_type.upper() != "I" or section_type.upper() != "U":
        section_type = "I"
    
    bf_top = float(input("Enter the b for top flange: "))
    tf_top = float(input("Enter the t for top flange: "))

    bf_bot = float(input("Enter the b for bottom flange: "))
    tf_bot = float(input("Enter the t for bottom flange: "))
        
    hw = float(input("Enter the hw: ")) 
    tw = float(input("Enter the t for web "))
    
    return Section(section_type.upper(), E,Fy,bf_top,bf_bot,tf_top,tf_bot,hw,tw)

def main():
    temp_section_props=""  
    temp_reaction =""
    temp_shear =""  
    temp_moment = ""
    temp_slope = ""
    temp_deflection = ""
    temp_check =""


    E = float(input("Enter E in mpa: "))
    Fy = float(input("Entet the fy in mpa: "))
    
    section = sectioncalc(E,Fy)
    temp_section_props =  section.print()
    
    n_spans = int(input("Enter the number of spans: "))
    
    spans_length = []
    for i in range(n_spans):
        L = float(input("Enter the {}th span length : ".format(i+1)))    
        spans_length.append(L)
    
    total_length = sum(spans_length)

    reactions_x = [0 for i in range(n_spans+1)]
    for i in range(1,len(reactions_x)):
        reactions_x[i] = sum(spans_length[:i])
    

    n_loads = int(input("Enter the number of loads: "))
    dist_loads = []
    pointed_loads = []
    total_load = 0

    for  i in range(n_loads):
        load_type = input("enter the type of load (1 for pointed 2 for distrbuted): ")
        try : 
            load_type = int(load_type)
        except:
            load_type = 1
        
        if load_type == 1 :
            x = float(input("Enter the load location from left : "))
            p = float(input("Enter the load value (+ for upward - for downward) : "))
            pointed_loads.append(SFunction(x,p,0))
            total_load += p
        else : 
            start = float(input("Enter the start location from left : "))
            end = float(input("Enter the end location from left : "))
            p = float(input("Enter the load value (+ for upward - for downward) : "))
            dist_loads.append(SFunction(start,p,0))
            dist_loads.append(SFunction(end,-1*p,0))
            total_load += (abs(end - start))*p
    


    reactions = []
    for x in reactions_x:
        reaction = Reaction(x)
        reactions.append(reaction)

    sfs_in_deflection = []
    for dist in dist_loads:
        sfs_in_deflection.append(dist.integrate(4))
    for pointed in pointed_loads:
        sfs_in_deflection.append(pointed.integrate(3))
    
    sfs_in_moment = []
    for dist in dist_loads:
        sfs_in_moment.append(dist.integrate(2))
    for pointed in pointed_loads:
        sfs_in_moment.append(pointed.integrate(1))

    
    total_factors = []
    consts = []
    for i in range(1,len(reactions_x)):
        x = reactions_x[i]
        const = 0
        for i in sfs_in_deflection:
            const += i.valueAt(x)
        consts.append(-1*const)

        factors = []
        for r in reactions:
            factors.append(r.sf_in_deflection.valueAt(x))
        factors.append(x)
        
        total_factors.append(factors)
    
    factor = []
    for r in reactions:
        factor.append(r.sf_in_moment.valueAt(x))
    factor.append(0)
    const = 0
    for i in sfs_in_moment:
        const += i.valueAt(x)
    consts.append(-1*const)

    total_factors.append(factor)

    l = [1 for i in range(len(reactions))]
    l.append(0)
    total_factors.append(l)
    consts.append(total_load*-1)


    reactions_values = np.linalg.solve(np.array(total_factors),np.array(consts))

    for i in range(len(reactions)):
        reactions[i].value = reactions_values[i]


    t = ""
    for r in reactions:
        t += r.print()
    temp_reaction = t


    shears = []
    for i in dist_loads:
        shears.append(i.integrate())
    for i in pointed_loads:
        shears.append(i)
    for reaction in reactions:
        shears.append(SFunction(reaction.cord,reaction.value,0))
    
    shears.sort(key=lambda x:x.start)
    
    t= ""
    for s in shears:
        t += s.print()
    temp_shear = t

    moments = []
    for i in shears:
        moments.append(i.integrate())
    moments.sort(key=lambda x:x.start)

    t=""
    for s in moments:
       t+= s.print()
    temp_moment = t

    slopes = []
    for i in moments:
        slopes.append(i.integrate())
    slopes.sort(key=lambda x:x.start)
    t = ""
    for s in slopes:
       t+= s.print()
    temp_slope = t
    
    deflections = []
    for i in slopes:
        deflections.append(i.integrate())
    deflections.append(SFunction(0,reactions_values[-1],1))
    deflections.sort(key=lambda x:x.start)
    
    t= ""
    for s in deflections:
       t+= s.print()
    temp_deflection=t
    
    x = 0
    xth = []
    shear_points = []
    while x < total_length:
        shear_point = 0
        for shear in shears :
            shear_point += shear.valueAt(x)
        shear_points.append(shear_point)
        xth.append(x)
        x +=0.001
    shear_points.append(0)
    xth.append(total_length)
    
    moment_points =[]
    for x in xth:
        moment_point = 0
        for moment in moments :
            moment_point += moment.valueAt(x)
        
        moment_points.append(moment_point)
    
    deflection_points = []
    for x in xth:
        deflection_point = 0
        for deflection in deflections :
            deflection_point += deflection.valueAt(x)
        deflection_points.append((deflection_point*1000000)/(E*section.I))
    
    plt.figure(0,dpi=200)
    plt.title("Shear")
    plt.xlabel("x (M)")
    plt.ylabel("Shear (KN)")
    plt.grid(True)
    plt.plot(xth,shear_points,label="nominal")
    plt.Axes.axis
    plt.axhline(linewidth=1, color='black')
    plt.savefig("shear.png")
    
    plt.figure(1,dpi=200)
    plt.title("Moment")
    plt.xlabel("x (M)")
    plt.axhline(linewidth=1, color='black')
    plt.ylabel("Moment (KN)")
    plt.grid(True)
    plt.plot(xth,moment_points,label="nominal")
    plt.Axes.axis
    plt.savefig("moment.png")


    plt.figure(2,dpi=200)
    plt.title("deflections")
    plt.xlabel("x (M)")
    plt.ylabel("deflection (mm)")
    plt.axhline(linewidth=1, color='black')
    plt.grid(True)
    plt.plot(xth,deflection_points,label="nominal")
    plt.Axes.axis
    plt.savefig("deflections.png")

   



    # ------------------------- checking the section-------------------------------

    for i in range(n_spans):
        
        temp_check +="<br><br><hr><br><br>Checking span {} : ".format(i+1)
        
        x_start = reactions_x[i]
        x_end = reactions_x[i+1]
        
        min_def = 1000
        start_index = -1
        for x in xth:
            deference = abs(x-x_start)
            if deference< min_def:
                start_index = xth.index(x)
                min_def = abs(x-x_start)

        min_def = 1000
        end_index = -1
        for x in xth[start_index:]:
            deference = abs(x-x_end)
            if deference< min_def:
                end_index = xth.index(x)
                min_def = abs(x-x_end)
        
        VUs = shear_points[start_index:end_index]
        MUs = moment_points[start_index:end_index]
        positive_MU = max(MUs)
        negative_MU = min(MUs)
        VU_max = max(max(VUs),abs(min(VUs)))

        temp_check += "<br><br>-------------- Shear --------------"
        Aw = section.h * section.tw
        c = math.sqrt(5*E/Fy)
        if section.h/section.tw <=1.1*c:
            Cv = 1
        elif  section.h/section.tw <=1.37*c:
            Cv = 1.1*c/(section.h/section.tw)
        else:
            Cv = 1.51*5*E/((section.h/section.tw)**2*Fy)
        
        Vn = 0.6*Cv*Aw*Fy

        temp_check +="<br><br>Aw = {:0.2f}, Cv = {:0.2f}".format(Aw,Cv)

        if VU_max <= 0.9*Vn/1000 :
            t ="<br><br>section meets the needs<br>0.9*Vn = {:0.2f} and Vu = {:0.2f}".format(0.9*Vn/1000 , VU_max)
        else :
            t= "<br><br>section doesnt meet the needs<br>0.9*Vn = {:0.2f} and Vu = {:0.2f}".format(0.9*Vn/1000 , VU_max)
        temp_check += t

        if section.section_type == "I":
            # ---- taslim bal feshari -----
            temp_check += "<br><br>-------------- taslim bal feshari --------------"

            if positive_MU > abs(negative_MU) :
                Sxc = section.Sx_top
                Iyc = section.I_top
            else : 
                Sxc = section.Sx_bot
                Iyc = section.I_bot
            
            Mp = min(Fy*section.Z,1.6*Fy*Sxc)
            Myc = Sxc*Fy

            if Iyc/section.I > 0.23:
                if section.landa <= section.landa:
                    Rpc = Mp/Myc
                else:
                    Rpc = min(Mp/Myc , ((section.landa - section.landa_p_web)/(section.landa_r_web-section.landa_p_web))*(Mp/Myc)*(Mp/Myc-1))
            else:
                Rpc = 1
            Mn_1 = Rpc*Sxc*Fy
            temp_check +="<br><br>Sxc = {:0.2f}, Iyc = {:0.2f}, MP = {:0.2f}, Myc = {:0.2f}, Rpc = {:0.2f}, Mn_1 = {:0.2f}".format(Sxc,Iyc,Mp,Myc,Rpc,Mn_1)

            # -------- kamanesh mozei bal feshari
            temp_check += "<br><br>-------------- khamesh mozei bal feshari --------------"

            if positive_MU > abs(negative_MU) :
                bf = section.bf_top
                tf = section.tf_top
                Sxc = section.Sx_top
                Sxt = section.Sx_bot
            else : 
                bf = section.bf_bot
                tf = section.tf_bot
                Sxt = section.Sx_top
                Sxc = section.Sx_bot
            
            Myc = Sxc*Fy
            landa = bf/(2*tf)
            Kc = max(0.35,float(4/math.sqrt(landa)))
            Kc = min(Kc,0.76)
            landa_p = 0.38*math.sqrt(E/Fy)
            if Sxt/Sxc >= 0.7:
                Fl = 0.7*Fy
            else:
                Fl = max(0.5*Fy, Sxt/Sxc*Fy)
            landa_r = 0.95*math.sqrt((Kc*E)/Fl)

            if landa <= landa_p:
                Mn_2 = Mn_1
            else:
                Myc = Sxc*Fy
                Mn_2 = Rpc*Myc-(Rpc*Myc-Fl*Sxc)*((landa-landa_p)/(landa_r-landa_p))
            
            temp_check +="<br><br>bf = {:0.2f}, tf = {:0.2f}, Sxc = {:0.2f}, Sxt = {:0.2f}, Myc = {:0.2f}, landa = {:0.2f}, KC = {:0.2f}, Fl = {:0.2f}, MN_2 = {:0.2f}".format(bf,tf,Sxc,Sxt,Myc,landa,Kc,Fl,Mn_2)
            
            # ------- taslim bal kesheshi
            temp_check +="<br><br>-------------- taslim bal kesheshi --------------"

            if positive_MU > abs(negative_MU) :
                Sxc = section.Sx_top
                Sxt = section.Sx_bot
                Iyc = section.I_top
            else : 
                Sxc = section.Sx_bot
                Sxt = section.Sx_top
                Iyc = section.I_bot
            
            Myt = Sxt*Fy
            if Sxt >=Sxc :
                Mn_3 = Mn_1
            else: 
                if landa <= landa_p:
                    Rpt = Mp/Myt
                else:
                    Rpt = min(Mp/Myt - ((Mp/Myt - 1)*((landa-landa_p)/(landa_r-landa_p)),Mp/Myt))
                Mn_3 = Rpt*Myt

            temp_check += "<br><br>Sxc = {:0.2f}, Sxt = {:0.2f}, Iyc = {:0.2f}, Myt = {:0.2f}, MN_3 ={:0.2f}".format(Sxc,Sxt,Iyc,Myt,Mn_3)

            Mn = min(Mn_1,Mn_2,Mn_3)
            Mu = max(positive_MU , abs(negative_MU))
            temp_check += "<br><br>------------------------------------------------------------"
            temp_check += "<br><br>Mn_1 = {:0.2f}<br> Mn_2 = {:0.2f}<br> Mn_3 = {:0.2f}".format(Mn_1/1000000,Mn_2/1000000,Mn_3/1000000)

            if Mu <= 0.9*Mn/1000000 :
                t = "<br><br>section meets the needs<br> 0.9*Mn = {:0.2f} and Mu = {:0.2f}".format(0.9*Mn/1000000 , Mu)
            else :
                t ="<br><br>section doesnt meet the needs<br> 0.9*Mn = {:0.2f} and Mu = {:0.2f}".format(0.9*Mn/1000000 , Mu)
            temp_check += t
        else:
            Mu = max(positive_MU , abs(negative_MU))
            Mn = section.Z * Fy

            if Mu <= 0.9*Mn/1000000 :
                t = "<br><br>section meets the needs<br> 0.9*Mn = {:0.2f} and Mu = {:0.2f}".format(0.9*Mn/1000000 , Mu)
            else :
                t ="<br><br>section doesnt meet the needs<br> 0.9*Mn = {:0.2f} and Mu = {:0.2f}".format(0.9*Mn/1000000 , Mu)
            temp_check += t

            
    
    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader("template"))
    template = jinja_env.get_template("template.html")
    jinja_var = {
        "temp_section_props":temp_section_props, 
        "temp_reaction":temp_reaction,
        "temp_shear":temp_shear,
        "temp_moment":temp_moment,
        "temp_slope":temp_slope,
        "temp_deflection":temp_deflection,
        "temp_check":temp_check
    }
    html = template.render(jinja_var)
    with open("report.html","w") as f:
        f.write(html)
    weasyprint.HTML("report.html").write_pdf("report.pdf")

    print("Your report is ready and stored in report.pdf")

if __name__ == "__main__":
    main()
