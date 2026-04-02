import numpy as np
from scipy import stats

def search(keyName, keyValue, listOfDicts):
    return filter(lambda item: item[keyName] == keyValue, listOfDicts)


def getEditingSitesList():
    editingSitesList = ['A->C', 'A->G', 'A->U', 
                        'C->A', 'C->G', 'G->U', 
                        'G->A', 'G->C', 'G->U', 
                        'U->A', 'U->C', 'U->G']
    
    return editingSitesList


"""
def getBoxPlotData(valuesList):
    
    # List sorting
    
    data = {}
    
    if valuesList:
        valuesList.sort()
        vList = np.array(valuesList)
        
        q1 = np.percentile(vList, 25, interpolation='linear')
        q3 = np.percentile(vList, 75, interpolation='linear')
        iqr = q3 - q1
        median = np.median(vList)
        lowerFence = q1 - 1.5*iqr
        upperFence = q3 + 1.5*iqr
        minimum = min(valuesList)
        maximum = max(valuesList)
        
        
        median = np.median(vList)
        upper_quartile = np.percentile(vList, 75)
        lower_quartile = np.percentile(vList, 25)
        
        iqr = upper_quartile - lower_quartile
        upper_whisker = vList[vList<=upper_quartile+1.5*iqr].max()
        lower_whisker = vList[vList>=lower_quartile-1.5*iqr].min()        
        
        lowerOutlier = []
        upperOutlier = []
        
        for value in vList:
            
            if value < lower_whisker:
                lowerOutlier.append(value)
                
            if value > upper_whisker:
                upperOutlier.append(value)        
        
        data = {
            'q1': lower_quartile,
            'q3': upper_quartile,
            'iqr': iqr,
            'median': median,
            'lowerFence': lower_whisker,
            'upperFence': upper_whisker,
            'lowerOutlier': lowerOutlier,
            'upperOutlier': upperOutlier,
            'minimum': minimum,
            'maximum': maximum
        }
    
    return data
"""  


def getBoxPlotData(vL):
    
    # List sorting
    
    data = {}
    
    valuesList = vL
    
    if valuesList:
        valuesList.sort()
        vList = np.array(valuesList)
        
        """
        q1 = np.percentile(vList, 25)
        q3 = np.percentile(vList, 75)
        iqr = q3 - q1
        median = np.median(vList)
        lowerFence = q1 - 1.5*iqr
        upperFence = q3 + 1.5*iqr
        upperOutlier = [ item for item in vList if item <(q1+1.5*iqr) ]
        lowerOutlier = [ item for item in vList if item >(q1-1.5*iqr) ]
        """

        
        median = np.median(vList)
        q3 = np.percentile(vList, 75, interpolation='linear')
        q1 = np.percentile(vList, 25, interpolation='linear')
        
        iqr = q3 - q1
        upper_whisker = vList[vList<=q3+1.5*iqr].max()
        lower_whisker = vList[vList>=q1-1.5*iqr].min()
        lowerFence = q1 - 1.5*iqr
        upperFence = q3 + 1.5*iqr        
        upperOutlier = []
        lowerOutlier = []
    
        for value in vList:
            
            if value < lower_whisker:
                lowerOutlier.append(value)
                
            if value > upper_whisker:
                upperOutlier.append(value)        
        
        data = {
            'q1': q1,
            'q3': q3,
            'iqr': iqr,
            'median': median,
            'lowerFence': lowerFence,
            'upperFence': upperFence,
            'lowerOutlier': lowerOutlier,
            'upperOutlier': upperOutlier,
            'upper_whisker': upper_whisker,
            'lower_whisker': lower_whisker
        }
    
    return data


def getColorPicker():
    
    colors = [
        '#3366cc', '#006666', '#009999',
        '#0000cc', '#0099cc', '#00ccff',
        '#003366', '#6600cc', '#3366ff',
        '#00cc99', '#33ccff', '#6666ff',
        '#6600cc', '#00cc66', '#99ccff',
        '#9933ff', '#00cc00', '#99ffcc',
        '#009933', '#ff00ff', '#cc9900',
        '#669900', '#ff9999', '#cc3399',
        '#cccc00', '#ff9933', '#cc0066',
        '#cc6600', '#ff0000', '#993333'
    ]
    
    return colors


def getRGBAColorPicker():
    
    colors = [
        'rgba(51, 102, 204, 0.4)', 
        'rgba(0, 102, 102, 0.4)', 
        'rgba(0, 153, 153, 0.4)',
        'rgba(0, 0, 204, 0.4)', 
        'rgba(0, 153, 204, 0.4)', 
        'rgba(0, 204, 255, 0.4)',
        'rgba(0, 51, 102, 0.4)', 
        'rgba(102, 0, 204, 0.4)', 
        'rgba(51, 102, 255, 0.4)',
        'rgba(0, 204, 153, 0.4)', 
        'rgba(51, 204, 255, 0.4)', 
        'rgba(102, 102, 255, 0.4)',
        'rgba(102, 0, 204, 0.4)', 
        'rgba(0, 204, 102, 0.4)', 
        'rgba(153, 204, 255, 0.4)',
        'rgba(153, 51, 255, 0.4)', 
        'rgba(0, 204, 0, 0.4)', 
        'rgba(153, 255, 204, 0.4)',
        'rgba(0, 153, 51, 0.4)', 
        'rgba(255, 0, 255, 0.4)', 
        'rgba(204, 153, 0, 0.4)',
        'rgba(102, 153, 0, 0.4)', 
        'rgba(255, 153, 153, 0.4)', 
        'rgba(204, 51, 153, 0.4)',
        'rgba(204, 204, 0, 0.4)', 
        'rgba(255, 153, 51, 0.4)', 
        'rgba(204, 0, 102, 0.4)',
        'rgba(204, 102, 0, 0.4)', 
        'rgba(255, 0, 0, 0.4)', 
        'rgba(153, 51, 51, 0.4)'
    ]
    
    return colors


def getRepositoryVersion():
    
    repos = {
        'TCGA': 'GRCh37.p5'
    }
    
    return repos