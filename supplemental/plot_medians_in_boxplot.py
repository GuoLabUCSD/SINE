
def plot_medians_in_boxplot(ax, df=None, x=None, y=None, hue=None, medians=None, vertical_offset=None, 
    hue_left_right=None, verbose=False, fontcolor='w', round_by=2):
    '''Adds medians to seaborn boxplot. Right now, hue only works for 2.
    example vertical_offset: 
    example hue_left_right: ['left_hue_var','right_hue_var']
    '''
    print(f'hue_left_right: {hue_left_right}, hue:{hue}')
    if all(x is None for x in [df, x, y, vertical_offset, medians]):
        print('You must provide either df/x/y or vertical_offset/medians')
        raise Error
    
    if medians is None:
        medians = df.groupby([x])[y].median()
    if vertical_offset is None:
        vertical_offset = df[y].median() * 0.01

    if hue_left_right is not None and hue is not None:
        medians = df.groupby([x,hue])[y].median()
        for xtick_i, xtick in enumerate(ax.get_xticklabels()):
            
            ax.text(xtick_i-.2, medians.loc[(xtick.get_text(),hue_left_right[0])]+vertical_offset, round(medians.loc[(xtick.get_text(),hue_left_right[0])],round_by),  horizontalalignment='center',  size='small', color=fontcolor, weight='semibold')
            ax.text(xtick_i+.2, medians.loc[(xtick.get_text(),hue_left_right[1])]+vertical_offset, round(medians.loc[(xtick.get_text(),hue_left_right[1])],round_by), horizontalalignment='center', size='small', color=fontcolor, weight='semibold')

    else:

        for xtick_i, xtick in enumerate(ax.get_xticklabels()):
            xtick = xtick.get_text()
            ax.text(xtick_i, medians[xtick] + vertical_offset, round(medians[xtick],round_by), 
                    horizontalalignment='center',size='small',color=fontcolor,weight='semibold')
        
        