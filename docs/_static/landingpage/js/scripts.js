var mr_firstSectionHeight,
    mr_nav,
    mr_navOuterHeight,
    mr_navScrolled = false,
    mr_navFixed = false,
    mr_outOfSight = false,
    mr_floatingProjectSections,
    mr_scrollTop = 0;
    mr_transforms               = ["transform","msTransform","webkitTransform","mozTransform","oTransform"],
    mr_transformproperty        = getSupportedPropertyName(mr_transforms),
    mr_scrollElement            = window,
   
$(document).ready(function(){
	

	"use strict"; 
	
	// Disable default behaviour for links that have # only as href
	
	$('a').click(function() {
        if ($(this).attr('href') === '#') {
            return false;
        }
    });
    
	
	// Open responsive menu
	
	$('.mobile-toggle').click(function(){
		$('nav').toggleClass('open-menu');
	});
		
	// Append .background-image-holder <img>'s as CSS backgrounds
	
	$('.background-image-holder').each(function(){
		var imgSrc= $(this).children('img').attr('src');
		$(this).css('background', 'url("' + imgSrc + '")');
    	$(this).children('img').hide();
        $(this).css('background-position', '50% 50%');
	});
	
	$('.foreground-image-holder').each(function(){
		var imgSrc= $(this).children('img').attr('src');
		$(this).css('background', 'url("' + imgSrc + '")');
    	$(this).children('img').hide();
        $(this).css('background-position', '50% 50%');
	});
	
	// Initialize YT player for youtube backgrounds
	
	$('.ytplayer').each(function(){
		var vidID = $(this).attr('data-video-id');
		$(this).attr('data-property','{videoURL:\''+vidID+'\',containment:\'.video-holder\',autoPlay:true, mute:true, startAt:28,opacity:1}');
	});
	
	
	// Set the height of the nav container to avoid jank on relative nav
	
	if($('nav').hasClass('relative-nav')){
		$('.nav-container').css('height', $('nav').outerHeight());
	}
		
	// Trigger menu toggle on nav-2
	
	$('.nav-2 .menu-toggle').click(function(){
		$(this).toggleClass('form-cross');
		$(this).closest('nav').find('.menu').toggleClass('expand');
	});
	
	$('.nav-3 .menu-toggle').click(function(){
		$(this).toggleClass('form-cross');
		$(this).closest('nav').find('.fullscreen-container').toggleClass('expand');
		setTimeout(function(){
			$('.fullscreen-container .vertical-align').toggleClass('show-menu');
		}, 300);
	});
	
	// Disable default behaviour for dummy links
	
	$('a').click(function(){
		if($(this).attr('href') == '#'){
			return false;
		}
	});
	
	// Nav Search input trigger
	
	$('.search-bar i').click(function(){
		$(this).closest('.search-bar').find('input[type="text"]').focus();
		if($(this).closest('.search-bar').find('input[type="text"]').val()){
			$(this).closest('.search-bar').find('input[type="submit"]').trigger('click');
		}
	});
	
	// Expanding Lists
	
	$('.expanding-list .title').click(function(){
		$(this).closest('li').toggleClass('active');
	});
	
	// Quantity control for products
	
	$('.add-to-cart .less').click(function(){
		if($(this).closest('form').find('input[type="text"]').val() >= 1){
			var newVal = parseInt($(this).closest('form').find('input[type="text"]').val()) - 1;
			$(this).closest('form').find('input[type="text"]').val(newVal);
		}
	});
	
	$('.add-to-cart .more').click(function(){
		var newVal = parseInt($(this).closest('form').find('input[type="text"]').val()) + 1;
		$(this).closest('form').find('input[type="text"]').val(newVal);
	});
    
    
	$('.hover-background').each(function(){
		$(this).mousemove(function( event ) {
			if(!mr_parallax.getScrollingState()){
				$(this).find('.background-image-holder').css(mr_transformproperty, 'translate3d(' + -event.pageX /18 + 'px,' + -(event.pageY-(mr_getScrollPosition())) /50+ 'px, 0)');
			}
				$(this).find('.layer-1').css(mr_transformproperty, 'translate3d(' + -event.pageX /30 + 'px,' + -event.pageY /30+ 'px, 0px)');
				$(this).find('.layer-2').css(mr_transformproperty, 'translate3d(' + -event.pageX /20 + 'px,' + -event.pageY /20+ 'px, 0px)');
			
		});
	});
    
   
    $('.promo-1 .btn').mouseenter(function(){
    	$(this).closest('.promo-1').find('.promo-image-holder').css('transform', 'scale(1.05)');
    }).mouseleave(function(){
    	$(this).closest('.promo-1').find('.promo-image-holder').css('transform', 'scale(1)');
    });
    
    // Project hover effects
    
    $('.hover-content').mouseenter(function(){
    	$(this).closest('.project').find('.background-image-holder').addClass('zoom');
    }).mouseleave(function(){
    	$(this).closest('.project').find('.background-image-holder').removeClass('zoom');
    });
    
    // Project filters
    
    // Populate filters
    $('.contained-projects').each(function(){
        
        var filters = "";

        $(this).find('.project').each(function(){
            
            var filterTags = $(this).attr('data-filter').split(',');
            
            filterTags.forEach(function(tagName){
                if(filters.indexOf(tagName) == -1){
                    filters += '<li data-filter="'+tagName+'">'+capitaliseFirstLetter(tagName)+'</li>';       
                }
            });
            $(this).closest('.contained-projects')
                   .find('ul.filters').empty().append('<li data-filter="all" class="active">All</li>').append(filters);
        });
    });

    $('.filters li').click(function(){
    	var filter = $(this).attr('data-filter');
    	$(this).closest('.filters').find('li').removeClass('active');
    	$(this).addClass('active');
    	
    	$(this).closest('.contained-projects').find('.project').each(function(){
            var filters = $(this).data('filter');
            
            if(filters.indexOf(filter) == -1){
                $(this).addClass('inactive');
            }
            else{
                $(this).removeClass('inactive');
            }
        });
    	
    	if(filter == 'all'){
    		$(this).closest('.contained-projects').find('.project').removeClass('inactive');
    	}
    });
    
    // Iframe fade-in dividers
    
    $('.video-strip .pre-video i').click(function(){
    	$(this).closest('.pre-video').addClass('fade-off');
    	$(this).closest('.video-strip').find('.iframe-holder').addClass('show-iframe');
    	var that = $(this);
    	setTimeout(function(){
    		that.closest('.video-strip').find('.iframe-holder').addClass('fade-on');
    	},500);
    });
    
    $('.video-strip .close-frame').click(function(){
    	$(this).closest('.iframe-holder').removeClass('fade-on');
    	var that = $(this);
    	setTimeout(function(){
    		that.closest('.video-strip').find('.iframe-holder').removeClass('show-iframe');
    		that.closest('.video-strip').find('.pre-video').removeClass('fade-off');
    	},500);
    });
    
    // Sliders
    
    $('.hero-slider').flexslider({ directionNav: false });
    $('.testimonials').flexslider({ directionNav: false });
    $('.image-slider').flexslider({ directionNav: false });
    
   // Twitter Feed
       jQuery('.tweets-feed').each(function(index) {
           jQuery(this).attr('id', 'tweets-' + index);
       }).each(function(index) {
           
           var TweetConfig = {
               "id": jQuery('#tweets-' + index).attr('data-widget-id'),
               "domId": '',
               "maxTweets": jQuery('#tweets-' + index).attr('data-amount'),
               "enableLinks": true,
               "showUser": true,
               "showTime": true,
               "dateFunction": '',
               "showRetweet": false,
               "customCallback": handleTweets
           };
           function handleTweets(tweets) {
               var x = tweets.length;
               var n = 0;
               var element = document.getElementById('tweets-' + index);
               var html = '<ul class="slides">';
               while (n < x) {
                   html += '<li>' + tweets[n] + '</li>';
                   n++;
               }
               html += '</ul>';
               element.innerHTML = html;
               return html;
           }
           twitterFetcher.fetch(TweetConfig);
       });
    
    // Instagram Feed
    
    if($('.instafeed').length){
    	jQuery.fn.spectragram.accessData = {
			accessToken: '1406933036.fedaafa.feec3d50f5194ce5b705a1f11a107e0b',
			clientID: 'fedaafacf224447e8aef74872d3820a1'
		};	

        $('.instafeed').each(function() {
            var feedID = $(this).attr('data-user-name') + '_';
            $(this).children('ul').spectragram('getUserFeed', {
                query: feedID,
                max: 12
            });
        });
    }   
    
    // Countdown
	
	$('.countdown').each(function(){
		$(this).countdown({until: new Date($(this).attr('data-date'))});
	});
    
    // Project Planner Form
    
    $('.planner-option').click(function(){
    	$(this).toggleClass('active');
    	if($(this).children('input').is(':checked')){
    		$(this).children('input').prop('checked', false);
    	}else{
    		$(this).children('input').prop('checked', true);
    	}
    });
    
    $('.planner-option,.planner-radio').each(function(){
    	var value = $(this).find('span').text();
    	$(this).find('input').attr('value', value);
    });
    
    $('.mock-radio').click(function(){
    	$(this).closest('.radio-group').find('.mock-radio').removeClass('active');
    	$(this).toggleClass('active');
    	if($(this).closest('.planner-radio').children('input').is(':checked')){
    		$(this).closest('.planner-radio').children('input').prop('checked', false);
    	}else{
    		$(this).closest('.planner-radio').children('input').prop('checked', true);
    	}
    });


    
    // Contact form code

    $('form.form-email, form.form-newsletter, form.project-planner').submit(function(e) {

        // return false so form submits through jQuery rather than reloading page.
        if (e.preventDefault) e.preventDefault();
        else e.returnValue = false;

        var thisForm = $(this).closest('form.form-email, form.form-newsletter, form.project-planner'),
            error = 0,
            originalError = thisForm.attr('original-error'),
            loadingSpinner, iFrame, userEmail, userFullName, userFirstName, userLastName, successRedirect;

        // Mailchimp/Campaign Monitor Mail List Form Scripts
        iFrame = $(thisForm).find('iframe.mail-list-form');

        thisForm.find('.form-error, .form-success').remove();
        thisForm.append('<div class="form-error" style="display: none;">' + thisForm.attr('data-error') + '</div>');
        thisForm.append('<div class="form-success" style="display: none;">' + thisForm.attr('data-success') + '</div>');


        if ((iFrame.length) && (typeof iFrame.attr('srcdoc') !== "undefined") && (iFrame.attr('srcdoc') !== "")) {

            console.log('Mail list form signup detected.');
            userEmail = $(thisForm).find('.signup-email-field').val();
            userFullName = $(thisForm).find('.signup-name-field').val();
            if ($(thisForm).find('input.signup-first-name-field').length) {
                userFirstName = $(thisForm).find('input.signup-first-name-field').val();
            } else {
                userFirstName = $(thisForm).find('.signup-name-field').val();
            }
            userLastName = $(thisForm).find('.signup-last-name-field').val();

            // validateFields returns 1 on error;
            if (validateFields(thisForm) !== 1) {
                console.log('Mail list signup form validation passed.');
                console.log(userEmail);
                console.log(userLastName);
                console.log(userFirstName);
                console.log(userFullName);

                iFrame.contents().find('#mce-EMAIL, #fieldEmail').val(userEmail);
                iFrame.contents().find('#mce-LNAME, #fieldLastName').val(userLastName);
                iFrame.contents().find('#mce-FNAME, #fieldFirstName').val(userFirstName);
                iFrame.contents().find('#mce-NAME, #fieldName').val(userFullName);
                iFrame.contents().find('form').attr('target', '_blank').submit();
                successRedirect = thisForm.attr('success-redirect');
                // For some browsers, if empty `successRedirect` is undefined; for others,
                // `successRedirect` is false.  Check for both.
                if (typeof successRedirect !== typeof undefined && successRedirect !== false && successRedirect !== "") {
                    window.location = successRedirect;
                }
            } else {
                thisForm.find('.form-error').fadeIn(1000);
                setTimeout(function() {
                    thisForm.find('.form-error').fadeOut(500);
                }, 5000);
            }
        } else {
            console.log('Send email form detected.');
            if (typeof originalError !== typeof undefined && originalError !== false) {
                thisForm.find('.form-error').text(originalError);
            }


            error = validateFields(thisForm);


            if (error === 1) {
                $(this).closest('form').find('.form-error').fadeIn(200);
                setTimeout(function() {
                    $(thisForm).find('.form-error').fadeOut(500);
                }, 3000);
            } else {
                // Hide the error if one was shown
                $(this).closest('form').find('.form-error').fadeOut(200);
                // Create a new loading spinner while hiding the submit button.
                loadingSpinner = jQuery('<div />').addClass('form-loading').insertAfter($(thisForm).find('input[type="submit"]'));
                $(thisForm).find('input[type="submit"]').attr('disabled', 'disabled');
                jQuery.ajax({
                    type: "POST",
                    url: "mail/mail.php",
                    data: thisForm.serialize()+"&url="+window.location.href,
                    success: function(response) {
                        // Swiftmailer always sends back a number representing numner of emails sent.
                        // If this is numeric (not Swift Mailer error text) AND greater than 0 then show success message.
                        $(thisForm).find('.form-loading').remove();

                        $(thisForm).find('input[type="submit"]').removeAttr('disabled');
                        if ($.isNumeric(response)) {
                            if (parseInt(response) > 0) {
                                // For some browsers, if empty 'successRedirect' is undefined; for others,
                                // 'successRedirect' is false.  Check for both.
                                successRedirect = thisForm.attr('success-redirect');
                                if (typeof successRedirect !== typeof undefined && successRedirect !== false && successRedirect !== "") {
                                    window.location = successRedirect;
                                }
                                thisForm.find('input[type="text"]').val("");
                                thisForm.find('textarea').val("");
                                thisForm.find('.form-success').fadeIn(1000);

                                thisForm.find('.form-error').fadeOut(1000);
                                setTimeout(function() {
                                    thisForm.find('.form-success').fadeOut(500);
                                }, 5000);
                            }
                        }
                        // If error text was returned, put the text in the .form-error div and show it.
                        else {
                            // Keep the current error text in a data attribute on the form
                            thisForm.find('.form-error').attr('original-error', thisForm.find('.form-error').text());
                            // Show the error with the returned error text.
                            thisForm.find('.form-error').text(response).fadeIn(1000);
                            thisForm.find('.form-success').fadeOut(1000);
                        }
                    },
                    error: function(errorObject, errorText, errorHTTP) {
                        // Keep the current error text in a data attribute on the form
                        thisForm.find('.form-error').attr('original-error', thisForm.find('.form-error').text());
                        // Show the error with the returned error text.
                        thisForm.find('.form-error').text(errorHTTP).fadeIn(1000);
                        thisForm.find('.form-success').fadeOut(1000);
                        $(thisForm).find('.form-loading').remove();
                        $(thisForm).find('input[type="submit"]').removeAttr('disabled');
                    }
                });
            }
        }
        return false;
    });

    $('.validate-required, .validate-email').on('blur change', function() {
        validateFields($(this).closest('form'));
    });

    $('form').each(function() {
        if ($(this).find('.form-error').length) {
            $(this).attr('original-error', $(this).find('.form-error').text());
        }
    });

    function validateFields(form) {
            var name, error, originalErrorMessage;

            $(form).find('.validate-required[type="checkbox"]').each(function() {
                if (!$('[name="' + $(this).attr('name') + '"]:checked').length) {
                    error = 1;
                    name = $(this).attr('name').replace('[]', '');
                    form.find('.form-error').text('Please tick at least one ' + name + ' box.');
                }
            });

            $(form).find('.validate-required').each(function() {
                if ($(this).val() === '') {
                    $(this).addClass('field-error');
                    error = 1;
                } else {
                    $(this).removeClass('field-error');
                }
            });

            $(form).find('.validate-email').each(function() {
                if (!(/(.+)@(.+){2,}\.(.+){2,}/.test($(this).val()))) {
                    $(this).addClass('field-error');
                    error = 1;
                } else {
                    $(this).removeClass('field-error');
                }
            });

            if (!form.find('.field-error').length) {
                form.find('.form-error').fadeOut(1000);
            }

            return error;
        }
    // End contact form code


	// Make map draggable when click
	
	$('.map-holder').click(function(){
		$(this).addClass('disable-overlay');
	});
	
	if($('.map-holder').length){
		$(window).scroll(function(){
			$('.map-holder').removeClass('disable-overlay');
		});
	}

	// Disable parallax on mobile

    if ((/Android|iPhone|iPad|iPod|BlackBerry|Windows Phone/i).test(navigator.userAgent || navigator.vendor || window.opera)) {
        $('section').removeClass('parallax');
    }

    // Load Google MAP API JS with callback to initialise when fully loaded
    if(document.querySelector('[data-maps-api-key]') && !document.querySelector('.gMapsAPI')){
        if($('[data-maps-api-key]').length){
            var script = document.createElement('script');
            var apiKey = $('[data-maps-api-key]:first').attr('data-maps-api-key');
            script.type = 'text/javascript';
            script.src = 'https://maps.googleapis.com/maps/api/js?key='+apiKey+'&callback=initializeMaps';
            script.className = 'gMapsAPI';
            document.body.appendChild(script);  
        } 
    }
   
    
});

///////////// End of Document Ready Function //////////////


$(window).load(function(){

	mr_fixedHeaderHeight        = $('.fixed-header').outerHeight(true);

	// Smooth scroll to inner links
	
	if($('nav').hasClass('nav-2')){
		$('.inner-link').smoothScroll({
			offset: -55,
			speed: 800
		});
	}else{
		var navHeight = $('nav').outerHeight();
		$('.inner-link').smoothScroll({
			offset: -navHeight,
			speed: 800
		});
    }
	
	// Initialize twitter feed
	
	var setUpTweets = setInterval(function(){
		if($('.tweets-slider').find('li.flex-active-slide').length){
			clearInterval(setUpTweets);
			return;
		}else{
			if($('.tweets-slider').length){
				$('.tweets-slider').flexslider({
					directionNav: false,
					controlNav: false
				});
			}
		}
    },500);
    
    // Show Background Images for sliders and dividers
    
    $('.background-image-holder').addClass('fadeIn');
    $('.foreground-image-holder').addClass('fadeIn');

    // Append Instagram BGs
    
    var setUpInstagram = setInterval(function(){
    	if($('.instafeed li').hasClass('bg-added')){
    		clearInterval(setUpInstagram);
			return;	
    	}else{
    		$('.instafeed li').each(function() {

				// Append background-image <img>'s as li item CSS background for better responsive performance
				var imgSrc = $(this).find('img').attr('src');
				$(this).css('background', 'url("' + imgSrc + '")');
				$(this).find('img').css('opacity', 0);
				$(this).css('background-position', '50% 0%');
				// Check if the slider has a color scheme attached, if so, apply it to the slider nav
				$(this).addClass('bg-added');
			});
			$('.instafeed').addClass('fadeIn');
    	}
    },500);
	
	
});

function capitaliseFirstLetter(string)
{
    return string.charAt(0).toUpperCase() + string.slice(1);
}


function getSupportedPropertyName(properties) {
    for (var i = 0; i < properties.length; i++) {
        if (typeof document.body.style[properties[i]] != "undefined") {
            return properties[i];
        }
    }
    return null;
}

function mr_getScrollPosition() {
        return mr_scrollElement != window ? mr_scrollElement.scrollTop : (document.documentElement.scrollTop == 0 ? document.body.scrollTop : document.documentElement.scrollTop);
    }


window.initializeMaps = function(){
    if(typeof google !== "undefined"){
        if(typeof google.maps !== "undefined"){
            $('.map-canvas[data-maps-api-key]').each(function(){
                    var mapInstance   = this,
                        mapJSON       = typeof $(this).attr('data-map-style') !== "undefined" ? $(this).attr('data-map-style'): false,
                        mapStyle      = JSON.parse(mapJSON) || [{"featureType":"landscape","stylers":[{"saturation":-100},{"lightness":65},{"visibility":"on"}]},{"featureType":"poi","stylers":[{"saturation":-100},{"lightness":51},{"visibility":"simplified"}]},{"featureType":"road.highway","stylers":[{"saturation":-100},{"visibility":"simplified"}]},{"featureType":"road.arterial","stylers":[{"saturation":-100},{"lightness":30},{"visibility":"on"}]},{"featureType":"road.local","stylers":[{"saturation":-100},{"lightness":40},{"visibility":"on"}]},{"featureType":"transit","stylers":[{"saturation":-100},{"visibility":"simplified"}]},{"featureType":"administrative.province","stylers":[{"visibility":"off"}]},{"featureType":"water","elementType":"labels","stylers":[{"visibility":"on"},{"lightness":-25},{"saturation":-100}]},{"featureType":"water","elementType":"geometry","stylers":[{"hue":"#ffff00"},{"lightness":-25},{"saturation":-97}]}],
                        zoomLevel     = (typeof $(this).attr('data-map-zoom') !== "undefined" && $(this).attr('data-map-zoom') !== "") ? $(this).attr('data-map-zoom') * 1: 17,
                        latlong       = typeof $(this).attr('data-latlong') != "undefined" ? $(this).attr('data-latlong') : false,
                        latitude      = latlong ? 1 *latlong.substr(0, latlong.indexOf(',')) : false,
                        longitude     = latlong ? 1 * latlong.substr(latlong.indexOf(",") + 1) : false,
                        geocoder      = new google.maps.Geocoder(),
                        address       = typeof $(this).attr('data-address') !== "undefined" ? $(this).attr('data-address').split(';'): false,
                        markerTitle   = "We Are Here",
                        isDraggable = $(document).width() > 766 ? true : false,
                        map, marker, markerImage,
                        mapOptions = {
                            draggable: isDraggable,
                            scrollwheel: false,
                            zoom: zoomLevel,
                            disableDefaultUI: true,
                            styles: mapStyle
                        };

                    if($(this).attr('data-marker-title') != undefined && $(this).attr('data-marker-title') != "" )
                    {
                        markerTitle = $(this).attr('data-marker-title');
                    }

                    if(address != undefined && address[0] != ""){
                            geocoder.geocode( { 'address': address[0].replace('[nomarker]','')}, function(results, status) {
                                if (status == google.maps.GeocoderStatus.OK) {
                                var map = new google.maps.Map(mapInstance, mapOptions); 
                                map.setCenter(results[0].geometry.location);
                                
                                address.forEach(function(address){
                                    var markerGeoCoder;
                                    
                                    markerImage = {url: window.mr_variant == undefined ? 'img/mapmarker.png' : '../img/mapmarker.png', size: new google.maps.Size(50,50), scaledSize: new google.maps.Size(50,50)};
                                    if(/(\-?\d+(\.\d+)?),\s*(\-?\d+(\.\d+)?)/.test(address) ){
                                        var latlong = address.split(','),
                                        marker = new google.maps.Marker({
                                                        position: { lat: 1*latlong[0], lng: 1*latlong[1] },
                                                        map: map,
                                                        icon: markerImage,
                                                        title: markerTitle,
                                                        optimised: false
                                                    });
                                    }
                                    else if(address.indexOf('[nomarker]') < 0){
                                        markerGeoCoder = new google.maps.Geocoder();
                                        markerGeoCoder.geocode( { 'address': address.replace('[nomarker]','')}, function(results, status) {
                                            if (status == google.maps.GeocoderStatus.OK) {
                                              
                                                marker = new google.maps.Marker({
                                                    map: map,
                                                    icon: markerImage,
                                                    title: markerTitle,
                                                    position: results[0].geometry.location,
                                                    optimised: false
                                                });
                                            }
                                        });
                                    }

                                });
                            } else {
                                console.log('There was a problem geocoding the address.');
                            }
                        });
                    }
                    else if(latitude != undefined && latitude != "" && latitude != false && longitude != undefined && longitude != "" && longitude != false ){
                        mapOptions.center   = { lat: latitude, lng: longitude};
                        map = new google.maps.Map(mapInstance, mapOptions); 
                        marker              = new google.maps.Marker({
                                                    position: { lat: latitude, lng: longitude },
                                                    map: map,
                                                    icon: markerImage,
                                                    title: markerTitle
                                                });

                    }

                }); 
        }
    }
}
initializeMaps();

// End of Maps


