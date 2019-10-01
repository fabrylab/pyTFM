def sum_on_area(masks,frame,res_dict,parameter_dict, db,db_info,label,x=None,y=None,sumtype="abs"):
    # check if mask it self is provided
    # masks is either a list of strings for clickpoints mask names or a mask as an array
    if isinstance(masks,np.ndarray):
        masks=masks.astype(bool)
        if sumtype=="abs":
            res_dict[frame]["%s "%(label)] = np.sum(np.sqrt(x[masks] ** 2 + y[masks] ** 2))
        if sumtype=="mean":
            res_dict[frame]["%s " % (label)] = np.mean(x[masks])
        if sumtype=="area": # area of original mask, without interpolation
            res_dict[frame]["%s " % (label)] = np.sum(masks)

    else:
        mtypes=make_iterable(masks)
        for mtype in mtypes:
            mask_membrane, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], mtype=mtype,
                                                   warn_thresh=1500)
            mask_membrane=cut_mask_from_edge(mask_membrane,parameter_dict["edge_padding"])
            if sumtype=="abs":
                mask_int = interpolation(binary_fill_holes(mask_membrane), dims=x.shape, min_cell_size=100)
                res_dict[frame]["%s on %s"%(label,mtype)] = np.sum(np.sqrt(x[mask_int] ** 2 + y[mask_int] ** 2))
            if sumtype=="mean":
                mask_int = interpolation(binary_fill_holes(mask_membrane), dims=x.shape, min_cell_size=100)
                res_dict[frame]["%s on %s" % (label, mtype)] = np.mean(x[mask_int])
            if sumtype=="area": # area of original mask, without interpolation
                area = np.sum(binary_fill_holes(mask_membrane)) * ((parameter_dict["pixelsize"] * 10 ** -6) ** 2)
                res_dict[frame]["%s of %s" % (label, mtype)] = area




def fill_patches_for_cell_layer(frame, parameter_dict,res_dict, db,db_info=None,single=True,**kwargs):
    # trying to load the mask from clickpoints
    print(db_info["frames_ref_dict"][frame])
    try:
        image=db.getImage(frame=db_info["frames_ref_dict"][frame]) #this is a work around, dont know why i cant use getMask directly
        mask = db.getMask(image).data
    # raising error if no mask object in clickpoints exist
    except AttributeError:
            raise Mask_Error("no mask of the cell membrane found for frame " + str(frame))
    mask_part = mask == 1  # type of "cell type1"


    mask_part = binary_fill_holes(mask_part).astype(bool)



    outer_lines=np.unique(labels[])
    mask[mask_part]=1
    mask[~mask_part] = 2
    mask=mask.astype(np.uint8)
    db.setMask(image=image,data=mask) # udapting mask in clickpoints
