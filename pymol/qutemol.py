import pymol


def qutemol(arg):
    '''

    DESCRIPTION

        "qutemol" QuteMol like image rendering settings

    USAGE
    
        qutemol on/off/0/1

        (c) Sergio M. Santos, University of Aveiro, 31st Jan. 2012
        
    '''
    settings = {'light_count'  :   10,
                'spec_count'   :    1,
                'shininess'    :   10,
                'specular'     : 0.25,
                'ambient'      :    0,
                'direct'       :    0,
                'reflect'      :  1.5,
                'depth_cue'    :    0,
                'bg_rgb'       : [1.0, 1.0, 1.0],
                'ray_shadow_decay_factor' :  0.1,
                'ray_shadow_decay_range'  :    2}
    
    if arg in ['on','On','ON','oN','1']:
        pymol.old_settings = {}
        for key in settings.keys():
            # backup previous settings
            pymol.old_settings[key] = cmd.get(key)
            # set new settings
            cmd.set(key, settings[key])
    else:
        # revert to previous settings
        for key in settings.keys():
            cmd.set(key, pymol.old_settings[key])
    return


cmd.extend('qutemol', qutemol)



def cartoonify(state='on'):
    '''

    DESCRIPTION

        "cartoonify" Cartoon like image rendering settings

    USAGE
    
        cartoonify on/off/0/1

        (c) Sergio M. Santos, University of Aveiro, 26th April 2012
        
    '''
    settings = {'light_count'   :     0,
                'ambient'       :     1,
                'line_radius'   : 0.006,
                'sphere_scale'  :  0.15,
                'ray_trace_mode':     1,
                'ray_trace_gain':     0,
                'bg_rgb'        : [1.0, 1.0, 1.0]}
    
    if state in ['on','On','ON','oN','1']:
        pymol.old_settings = {}
        for key in settings.keys():
            # backup previous settings
            pymol.old_settings[key] = cmd.get(key)
            # set new settings
            cmd.set(key, settings[key])
    else:
        # revert to previous settings
        for key in settings.keys():
            cmd.set(key, pymol.old_settings[key])
    return


cmd.extend('cartoonify', cartoonify)

