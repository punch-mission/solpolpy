class Depolarize:
    '''Depolarize - Transforms between different Polarization states
    
    The Depolarize module defines a depolarize object that represents an 
    transformation between different polarization states.  Depolarizer objects 
    can be used to transform between different polarization states.  The base 
    package is supplied with subclasses that implement several general-purpose 
    polarization transformations; additional subpackages supply suites of 
    polarization for different instruments.
    
    The simplest way to use a Depolarizer object is to transform polarized data
    frames between polarization representations.  The "apply" method accepts 
    multiple FITS, an ndcube array or dictionary of angles and dataframes 
    and transforms the dataframes according to the formulae embedded in the 
    object. 
    
    The output is an NDarray of dataframes representing the transformed 
    polarization frames.
    
    NOTE: The number of input data frames define which output polarization 
    representations can be generated. 3 polarization angles are required
    to generate an M, Z, P three-polarizer measurement and represenation 
    system. Four polarization angles are required to acurately derive the 
    full set of Stokes parameters.
    
    NOTE: Depolarizer considers images to be 2-D arrays indexed in conventional, 
    sane order: the pixel coordinate system is defined so that (0,0) is at the 
    *center* of the LOWER, LEFT pixel of an image, with (1,0) being one pixel 
    to the RIGHT and (0,1) being one pixel ABOVE the origin, i.e. pixel vectors
    are considered to be (X,Y) by default. This indexing method agrees with 
    nearly the entire scientific world aside from the NumPy community, which
    indexes image arrays with (Y,X), and the SciPy community, which sometimes
    indexes images arrays with (Y,-X) to preserve handedness.  For this reason,
    if you use Transformed vectors directly to index array data (outside of
    the map method) then you must reverse the order of normal Transform vector
    components.  A handy ArrayIndex subclassed Transform is supplied, to do 
    this conversion from sane coordinates to NumPy array index coordinates.
     
    Examples
    --------
        import depolarizer as d                          # Load the package
        PolarizationObject = d.depolarize(dictionary)    # Generate a depolarization object
        newData = PolarizationObject.apply(BpB)          # Apply the scaling to the object

    
    Methods
    -------
    
    Depolarization objects have the following public methods:
        
        - apply: accepts a polarizarion state to transform the data into.
                             
        - compose: returns a depolarization state that is the composition of this
            one with one or more existing depolarization state, 
            i.e. B, pB (input) -> M,Z,P -> Stokes
        
         
    Built-in Subclasses
    -------------------
    
    The Depolarize class is a container for the subclasses that do the actual 
    mathematical work.  Each subclass represents a depolarization operation; 
    an instance of that subclass represents a single depolarization selected 
    from the family by the depolarizers you pass to the constructor.
    
    The depolarizer module itself defines the following subclasses:
                    
        - Composition: groups together multiple depolarization transforms 
        into one composite depolarization transform.
          
        
    '''
    #def __init__(self):
    #   raise AssertionError(\
    #       "generic depolarizations must be subclassed (e.g. Depolarize.identity)"\
    #           )

    def __str__(self):

        '''
        __str__ - stringify a generic transform
        
        The Depolarize stringifier handles putting any generic descriptors
        on the output string. It wraps a more specific string that is 
        supplied by subclassed __str__s.  Because __str__ can't take any
        arguments directly, subclasses should store their output in 
        self._strtmp, then call super().__str__().  The ._strtmp gets
        deletedh ere.
        
        Some subclasses want to generate strings for additional Depolarizations 
        embedded inside them. Stringifying those sub-objects should be more 
        terse than a regular stringification, so there's a separate 
        ._str_not_top_tmp flag that they set to turn off the contextual portion 
        of the string.
        
        Returns
        -------
        str
            the string.
        '''
        try:
            s=self._strtmp
        except:
            return "Generic Depolarize stringifier - you should never see this"
        del self._strtmp

        # Some subclasses stringify by listing subsidiary Depolarizations. 
        # They call __str__ on those depolarizations, but want the verbose
        # "Depolarize()" suppressed. Those subclasses can set the flag in 
        # their subsidiary objects. If the flag is set then we return just 
        # the string we're passed. If the flag exists we always end up 
        # deleting it.
        try:
            flag = self._str_not_top_tmp
            del self._str_not_top_tmp
            if(flag):
                return s
        except:
            pass
        
        return f"Depolarize( {s} )"

    def resolve(self, out_polarize_state, separation=None, alpha=None, Error=False):
       
        '''
        Apply - apply a depolarization transform to a set of input
        dataframes.

        Parameters
        ----------
        out_polarize_state : string
          This is the polarization state you want to convert your input 
          dataframes to.  These include:
        
        - Stokes: convert the input dataframes to output Stokes dataframes. If
            3 polarization angles are input, or can be derived, the I, Q, and U
            Stokes parameters are output. If 4 or more polarization angles are
            provided, and the angles are conducsive, the full I,Q,U, and V 
            Stokes parameters are provided.

        - B: converts input dataframes into B ("unpolarized 
            brightness") parameters.

        - pB: converts input dataframes into pB ("polarized brightness") 
            parameters.

        - Bt: Produces "Bt", the radiance observed through a linear polarizer 
            oriented tangentially to a solar-concentric circle passing through 
            the image point of interest.

        - Br: Produces "Br" the radiance observed through a linear polarizer 
            oriented radially to the centre Sun through the same point.
                
        - BrBt: Produces both the  radiance observed through a linear polarizer 
            oriented tangentially to a solar-concentric circle passing through 
            the image point of interest "Bt", and the radiance observed through
            a linear polarizer oriented radially to the centre Sun through the 
            same point "Br".
            
        - 3pol or MZP: converts input dataframes into a system of virtual 
            polarizer triplets each separated by 60 degrees (Minus [-60], 
            Zero [0], Plus [60]).

        - 4pol: converts input dataframes into a system of virtual polarizer 
            triplets each separated by 45 degrees at -45, 0, 45, 90.
          
        - Xpol: converts input dataframes into a system of virtual polarizer 
            triplets each separated by X degrees, specified by the optional 
            input separation, which should be input in degrees.

        - BpB: converts input dataframes into two dataframes, B ("unpolarized 
            brightness") parameter and its counterpart pB ("polarized brightness")
        
        Raises
        ------
        AssertionError
          This gets raised if the data cannot be converted or polarization 
          transformation cannot calculated due to a discontinuity or infinity.
          
        Returns
        -------
        numpy.ndarray
            The transformed vector data are returned as a numpy.ndarray.  Most 
            Transforms maintain the dimensionality of the source vectors.  Some 
            embed (increase dimensionality of the vectors) or project (decrease
            dimensionality of the vectors); additional input dimensions, if 
            present, are still appended to the output vectors in all any case. 
        '''


    def depolarize(self, data):
        '''
        depolarize - input the data to be depolarized or converted
        
        Parameters
        ----------
        data : ndarray
          This is the data to which the Depolarization should be applied.
          The data should be in one of the following formats:
          
        - an NDcube object

        - a dictionary of 3 or more 
          
           NumPy ndarray, with the final axis 
          running across vector dimension.  Earlier dimensions are broadcast
          (e.g., a WxHx2 NumPy array is treated as a WxH array of 2-vectors).
          The -1 axis must have sufficient size for the transform to work.
          If it is larger, then subsequent vector dimensions are ignored , so 
          that (for example) a 2-D Transform can be applied to a WxHx3 NumPy 
          array and the final WxH plane is transmitted unchanged.  
          
        invert : Boolean, (default False)
          This is an optional flag indicating that the inverse of the transform
          is to be applied, rather than the transform itself. 
          
        Raises
        ------
        ValueError
          Dimensional mismatch or non-array inputs cause this to be thrown.
          
        AssertionError
          This gets raised if the Transform won't work in the intended direction.
          That can happen because some mathematical transforms don't have inverses;
          Transform objects contain a flag indicating validity of forward and 
          reverse application.
          
        Returns
        -------
        numpy.ndarray
            The transformed vector data are returned as a numpy.ndarray.  Most 
            Transforms maintain the dimensionality of the source vectors.  Some 
            embed (increase dimensionality of the vectors) or project (decrease
            dimensionality of the vectors); additional input dimensions, if 
            present, are still appended to the output vectors in all any case.
        '''
        # Start by making sure we have a dictionary
        #if( not isinstance(data, np.ndarray) ):
        #    data = np.array(data)
        data = self.apply(data)
            
        return data