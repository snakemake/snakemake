'use strict';

class MetaData extends React.Component {

    render() {
        let metadatalist = []
        // If workflow description is given, show it on the landing page with the metadata
        if (workflow_desc) {
            metadatalist.push(...
            [
                e(
                    "div",
                    {
                        className: "prose prose-sm max-w-lg",
                        dangerouslySetInnerHTML: { __html: workflow_desc }
                    }
                ),
            ]
            );
        }

        // Iterate over the metadata dictionary and add the entries to the landing page
        for (const [key, value] of Object.entries(metadata)) {
            metadatalist.push(... this.innerRender(key, value));
        }
        console.log(workflow_desc)
        return e(
            //workflow_desc,
            "ul",
            { },
            metadatalist
        )
    }

    innerRender(key, value) {
        // check if the value is an dictionary, if this is the case, recursively
        // add the elements to the landing page
        if (value.constructor == Object) {
            let result = [];
            for (const [inner_key, inner_value] of Object.entries(value)) {
                result.push(...this.innerRender(inner_key, inner_value));
            }
            return result;
        }
        return [
            e(
                ListHeading,
                { key: `${key}-heading`, text: key }
            ),
            e(
                ListItem,
                {
                    key: `${key}`,
                    className: "p-1",
                    // this is safe as the value is already rendered by docutils
                    dangerouslySetInnerHTML: { __html: value }
                }
            )
        ]
    }

}
